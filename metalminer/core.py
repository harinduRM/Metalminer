from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional
from collections import deque
import multiprocessing as mp
import os
import pickle

import numpy as np
import pandas as pd
import py3Dmol
import re

from ccdc import io, molecule, search
from scipy.spatial import cKDTree

from .auto_disorder import clean_small_cluster
from .constants import ELEMENT_NAMES
from .ligands import analyze_ligands
from .oxidation import (
    HAS_GEMINI,
    calculate_bvs,
    extract_oxidation_state_with_gemini,
    extract_text_oxidation_state,
    resolve_mixed_valence,
)
from .utils import (
    distance,
    generate_enriched_xyz_string,
    get_abstract_from_openalex,
    get_nuclearity_label,
)
from .validation import validate_geometry


@dataclass
class Config:
    TARGET_METAL_list: List[str] = field(default_factory=lambda: ["Pu"])
    target_csd_ids: Optional[List[str]] = None
    EXTRACTION_METHOD: str = "Topological"
    R_FACTOR_LIMIT: float = 10
    PROCESS_LIMIT: int = 10000
    SITE_TIMEOUT_SECONDS: Optional[int] = None
    RETRY_TIMEOUTS: bool = True
    VISUALIZE: bool = False
    VISUALIZATION_LIMIT: int = 15
    FILTER_POLYMERIC: bool = True
    FILTER_POWDER: bool = True
    FILTER_ALL_DISORDER: bool = False
    FILTER_PRIMARY_DISORDER: bool = False
    CORRECT_PRIMARY_DISORDER: bool = True
    Hydrogen_Addition_method: str = "Geometric"
    DISORDER_RESOLVE_METHOD: str = "Hybrid"
    num_metal_layers: int = 1
    OXIDATION_STATES_FILTER: List[str] = field(default_factory=lambda: ["all"])
    GET_ABSTRACT: bool = True
    Edit_manual: bool = True
    Geometric_radius: float = 3.8
    Extraction_Cut_off_distances: Dict[str, float] = field(
        default_factory=lambda: {"LIMIT_NM_NM": 2.8, "LIMIT_M_NM": 3.5, "LIMIT_H_X": 1.3}
    )
    metalloligands: List[str] = field(default_factory=lambda: ["Cr", "V", "Mo", "W", "S", "Co"])


@dataclass
class Result:
    results: List[Dict[str, object]]
    dataframe: Optional[pd.DataFrame]
    summary_report: str
    artifacts: Dict[str, Optional[str]]
    visualizations_html: List[str]


def default_emit(text: str) -> None:
    print(text)


class SimpleHit:
    """Helper class to make direct CSD entries mimic search hits."""

    def __init__(self, entry):
        self.entry = entry
        self.identifier = entry.identifier


RESULT_COLUMNS = (
    "CSD ID",
    "Site Label",
    "Metal",
    "Oxidation State",
    "Chemical Name",
    "Unit_cell",
    "Abstract",
    "Structure Type",
    "Coordination Number",
    "Ligand_count",
    "AN_n_Ligand_Info",
    "VALIDATION_FAILED",
    "HAD_RDKIT_ISSUE",
    "Bridging_Ligand",
    "Total_Ligand_Charge",
    "OS_Method",
    "Hydrogen Method",
    "OS_source",
    "AN_n_XYZ Coordinates",
    "Had_prismary_sphere_disorder",
    "DOI",
)
TIMEOUT_PREFIX = "Timeout::"


def _build_timeout_value(csd_id: str, site_label: str) -> str:
    return f"{TIMEOUT_PREFIX}{csd_id}::{site_label}"


def _parse_timeout_site(value: object) -> Optional[str]:
    if isinstance(value, str) and value.startswith(TIMEOUT_PREFIX):
        parts = value.split("::", 2)
        if len(parts) == 3:
            return f"{parts[1]}_{parts[2]}"
    if isinstance(value, str) and value.strip().lower() == "timeout":
        return "timeout"
    return None


def _process_site_for_timeout(payload):
    try:
        entry_reader = io.EntryReader("CSD")
        entry = entry_reader.entry(payload["csd_id"])
        if not entry:
            return {"status": "entry_missing", "error": "Entry not found in CSD."}

        config = payload["config"]
        target_metal = payload["target_metal"]
        clean_base_label = payload["clean_base_label"]

        core_result = process_structure_core(
            entry,
            target_metal,
            target_atom_label=clean_base_label,
            num_metal_layers=config.num_metal_layers,
            perform_correction=config.CORRECT_PRIMARY_DISORDER,
            EXTRACTION_METHOD=config.EXTRACTION_METHOD,
            disorder_resolve_method=config.DISORDER_RESOLVE_METHOD,
            Extraction_Cut_off_distances=config.Extraction_Cut_off_distances,
            Geometric_radius=config.Geometric_radius,
            metalloligands=config.metalloligands,
            emit=default_emit,
        )
        mol, center_metal, has_primary_disorder, nuclearity_count, bridging_indices, validation_failed = core_result

        if mol is None:
            if has_primary_disorder and config.CORRECT_PRIMARY_DISORDER:
                return {"status": "fail_correction"}
            return {"status": "fail_extraction"}

        if has_primary_disorder and config.FILTER_PRIMARY_DISORDER:
            return {"status": "skip_primary_disorder"}

        name_str = entry.chemical_name if entry.chemical_name else ""
        form_str = entry.formula if entry.formula else ""
        combined_text = f"{name_str} {form_str}"

        AN_n_ligand_data, has_bridging, o_map_n, virtual_hs = analyze_ligands(
            mol,
            center_metal,
            combined_text,
            bridging_indices,
            Hydrogen_Addition_method=config.Hydrogen_Addition_method,
        )

        is_valid = validate_geometry(
            mol,
            center_metal,
            virtual_hydrogens=virtual_hs,
            check_rdkit=True,
            emit=default_emit,
        )
        validation_failed = not is_valid

        AN_n_xyz_block = generate_enriched_xyz_string(mol, virtual_hs)

        bvs_val = calculate_bvs(center_metal, center_metal.neighbours)
        bvs_int = round(bvs_val)
        final_state = "?"
        source = "Unknown"
        source_content = ""
        method_used = ""

        candidate_states_list = []

        abstract_text = "N/A"
        should_fetch = config.GET_ABSTRACT or HAS_GEMINI

        if should_fetch and entry.publication and entry.publication.doi:
            doi_key = str(entry.publication.doi).strip().lower()
            fetched_abstract = get_abstract_from_openalex(doi_key)
            if fetched_abstract:
                abstract_text = fetched_abstract

        text_states, combined_text_os = extract_text_oxidation_state(entry, target_metal)
        if text_states:
            candidate_states_list.extend(text_states)
            source = "Text Mining"
            method_used = "Text Mining"
        elif HAS_GEMINI and abstract_text != "N/A":
            llm_states = extract_oxidation_state_with_gemini(abstract_text, target_metal)
            if llm_states:
                candidate_states_list.extend(llm_states)
                source = "Gemini LLM"
                method_used = "Gemini LLM"

        if len(candidate_states_list) > 1:
            final_state = resolve_mixed_valence(center_metal, candidate_states_list, mol)
            method_used += " (Resolved by Geometry)"
        elif len(candidate_states_list) == 1:
            final_state = f"+{candidate_states_list[0]}"
            source = source or "Text Mining"
            method_used = method_used or "Text Mining"
        else:
            final_state = f"+{bvs_int}"
            source = "BVS"
            method_used = "BVS (Geometry)"

        should_keep = False
        if not config.OXIDATION_STATES_FILTER or config.OXIDATION_STATES_FILTER == ["all"]:
            should_keep = True
        else:
            for allowed_state in config.OXIDATION_STATES_FILTER:
                if allowed_state.lower() == "actinyl":
                    if "actinyl" in (source_content or combined_text).lower():
                        should_keep = True
                        break
                else:
                    try:
                        allowed_int = int(str(allowed_state).replace("+", ""))
                        if final_state.lstrip("+") == str(allowed_int):
                            should_keep = True
                            break
                    except ValueError:
                        continue

        if not should_keep:
            return {"status": "skip_os_mismatch"}

        structure_type = get_nuclearity_label(nuclearity_count)
        coordination_number = sum(len(x["Connecting_atom"]) for x in AN_n_ligand_data)
        total_ligand_charge = sum(item["Charge"] * item["number"] for item in AN_n_ligand_data)
        lig_count = sum(x["number"] for x in AN_n_ligand_data)
        has_rdkit_issue = any(item.get("HAD_RDKIT_ISSUE", False) for item in AN_n_ligand_data)

        doi_str = entry.publication.doi if entry.publication else "N/A"

        uc_dict = {}
        if entry.crystal:
            uc_dict = {
                "a": entry.crystal.cell_lengths.a,
                "b": entry.crystal.cell_lengths.b,
                "c": entry.crystal.cell_lengths.c,
                "alpha": entry.crystal.cell_angles.alpha,
                "beta": entry.crystal.cell_angles.beta,
                "gamma": entry.crystal.cell_angles.gamma,
            }

        result_row = {
            "CSD ID": entry.identifier,
            "Site Label": clean_base_label,
            "Metal": target_metal,
            "Oxidation State": final_state,
            "Chemical Name": name_str,
            "Unit_cell": uc_dict,
            "Abstract": abstract_text,
            "Structure Type": structure_type,
            "Coordination Number": coordination_number,
            "Ligand_count": lig_count,
            "AN_n_Ligand_Info": AN_n_ligand_data,
            "VALIDATION_FAILED": validation_failed,
            "HAD_RDKIT_ISSUE": has_rdkit_issue,
            "Bridging_Ligand": has_bridging,
            "Total_Ligand_Charge": total_ligand_charge,
            "OS_Method": method_used,
            "Hydrogen Method": config.Hydrogen_Addition_method,
            "OS_source": source,
            "AN_n_XYZ Coordinates": AN_n_xyz_block,
            "Had_prismary_sphere_disorder": has_primary_disorder,
            "DOI": doi_str,
        }
        return {
            "status": "ok",
            "result_row": result_row,
            "method_used": method_used,
            "has_primary_disorder": has_primary_disorder,
        }
    except Exception as exc:
        return {"status": "error", "error": str(exc)}


def _process_site_worker(payload, output_queue):
    output_queue.put(_process_site_for_timeout(payload))


def _run_site_with_timeout(payload, timeout_seconds):
    ctx = mp.get_context("spawn")
    output_queue = ctx.Queue(maxsize=1)
    process = ctx.Process(target=_process_site_worker, args=(payload, output_queue))
    process.start()
    process.join(timeout_seconds)
    if process.is_alive():
        process.terminate()
        process.join()
        return {"status": "timeout"}
    if not output_queue.empty():
        return output_queue.get()
    return {"status": "error", "error": "Worker returned no result."}

def process_structure_core(entry, target_element, target_atom_label, num_metal_layers=1,
                            metalloligands=["Cr", "V", "Mo", "W", "S","Co"], 
                            perform_correction=True, EXTRACTION_METHOD="Topological",
                            disorder_resolve_method="Hybrid", Extraction_Cut_off_distances ={'LIMIT_NM_NM':2.8,'LIMIT_M_NM' :3.5,'LIMIT_H_X':1.3},
                            manual_commands=None, Geometric_radius = 3.8, emit: Callable[[str], None] = default_emit):
    """
    Extracts cluster around `target_atom_label` with strict connectivity filters.
    
    FIXES:
    - Enforces Hard Max Limits on learned bonds to reject "Long Contacts".
    - NonMetal-NonMetal limit: 2.8 A (rejects Cl...C 3.9 A).
    - Metal-Ligand limit: 3.5 A.
    - H-NonMetal limit: 1.3 A.
    """
    if not entry.has_3d_structure: return None, None, False, 0

    # --- 0. DEFINE CATEGORIES ---
    # Heuristic: If it's not a metal, treat as Non-Metal for bond limiting
    # (Note: CCDC 'is_metal' covers transition, lanthanides, actinides, etc.)

    excluded_metals = {"Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Sr", "Ba", "Al"}
    # Create a local set for fast lookup, removing the target if present
    active_metalloligands = set(metalloligands)
    if target_element in active_metalloligands:
        active_metalloligands.remove(target_element)
    
 
    
    # --- 1. LEARN BONDING (With Sanitation) ---
    ref_mol = entry.molecule
    valid_bond_ranges = {}
    
    default_cutoffs = {"LIMIT_NM_NM": 2.8, "LIMIT_M_NM": 3.5, "LIMIT_H_X": 1.3}
    if EXTRACTION_METHOD == "Topological" and Extraction_Cut_off_distances:
        cutoffs = Extraction_Cut_off_distances
    else:
        cutoffs = default_cutoffs

    # HARD LIMITS to reject artifacts
    LIMIT_NM_NM = cutoffs["LIMIT_NM_NM"]  # NonMetal-NonMetal Max (I-I is ~2.7)
    LIMIT_M_NM = cutoffs["LIMIT_M_NM"]  # Metal-Ligand Max (Actinide-Cl ~2.6-2.9)
    LIMIT_H_X = cutoffs["LIMIT_H_X"]  # Hydrogen-Any Max
    
    for bond in ref_mol.bonds:
        try:
            d = bond.length
            if d > 5.0: continue # Skip gross errors
            
            a1, a2 = bond.atoms[0], bond.atoms[1]
            s1, s2 = sorted([a1.atomic_symbol, a2.atomic_symbol])
            m1, m2 = a1.is_metal, a2.is_metal
            
            # --- FILTER 1: Reject Long Contacts Immediately ---
            if 'H' in (s1, s2):
                if d > LIMIT_H_X: continue
            elif not m1 and not m2:
                # Both Non-Metals (e.g. C-Cl, N-C, C-C)
                if d > LIMIT_NM_NM: continue
            elif m1 != m2:
                # Metal-Ligand (e.g. Np-Cl)
                if d > LIMIT_M_NM: continue
            
            # If passed, record it
            key = (s1, s2)
            if key not in valid_bond_ranges: valid_bond_ranges[key] = []
            valid_bond_ranges[key].append(d)
        except: continue

    bond_limits = {}
    tolerance = 0.20 
    
    global_max = 0.0
    for (s1, s2), lengths in valid_bond_ranges.items():
        max_d = max(lengths) + tolerance
        
        # --- FILTER 2: Clamp Learned Limits (Redundant Safety) ---
        # Even if a bad bond sneaked in, clamp the final limit used for injection
        if 'H' in (s1, s2):
            if max_d > LIMIT_H_X: max_d = LIMIT_H_X
        elif (s1 not in ELEMENT_NAMES and s2 not in ELEMENT_NAMES): # Approx non-metal check using your dict or is_metal logic
             # Note: Using strict 2.8 again here is safer
             if max_d > LIMIT_NM_NM: max_d = LIMIT_NM_NM
             
        bond_limits[(s1, s2)] = max_d
        bond_limits[(s2, s1)] = max_d
        if max_d > global_max: global_max = max_d
    
    # Clamp global max to prevent KDTree from searching too far
    if global_max < 4.5: global_max = 4.5
    if global_max > 6.0: global_max = 6.0 

    # --- 2. EXPAND ---
    try:
        cell_params = entry.crystal.cell_lengths
        min_dim = min(cell_params.a, cell_params.b, cell_params.c)
        est_radius = 12.0 + ((num_metal_layers - 1) * 10.0)
        n_cells = int(np.ceil(est_radius / min_dim))
        if num_metal_layers <=2: n_cells = 1 
        else: n_cells=2
        packed_mol = entry.crystal.packing(box_dimensions=((-n_cells, -n_cells, -n_cells), (n_cells, n_cells, n_cells)))
    except Exception as e:
        emit(f"Expansion failed for {entry.identifier}: {e}. Falling back to molecule.")
        packed_mol = entry.molecule

    # --- 3. INJECT CONNECTIVITY ---
    all_atoms = packed_mol.atoms
    coords = np.array([a.coordinates for a in all_atoms])
    symbols = np.array([a.atomic_symbol for a in all_atoms])
    
    tree = cKDTree(coords)
    
    unique_in_structure = np.unique(symbols)
    elem_to_int = {sym: i for i, sym in enumerate(unique_in_structure)}
    
    is_metal_map = {}
    for atom in all_atoms:
        if atom.atomic_symbol not in is_metal_map:
            is_metal_map[atom.atomic_symbol] = atom.is_metal

    n_types = len(unique_in_structure)
    cutoff_matrix = np.zeros((n_types, n_types))
    
    for i, sym1 in enumerate(unique_in_structure):
        for j, sym2 in enumerate(unique_in_structure):
            limit = bond_limits.get((sym1, sym2), 0.0)
            
            m1 = is_metal_map.get(sym1, False)
            m2 = is_metal_map.get(sym2, False)

            # Default logic for Metal-Ligand ONLY
            if limit == 0.0:
                # ONLY apply default if it is Metal-Ligand (m1 != m2)
                # AND strictly NOT H.
                if (m1 != m2) and (sym1 != 'H' and sym2 != 'H'):
                    limit = 3.5 # Safe default, reduced from 4.25
            
            # Force clamp on matrix entries to prevent any "leak"
            if not m1 and not m2 and limit > LIMIT_NM_NM:
                limit = LIMIT_NM_NM # Enforce 2.8 on NM-NM

            cutoff_matrix[i, j] = limit
            
    raw_pairs = tree.query_pairs(r=global_max)
    
    if raw_pairs:
        pairs_array = np.array(list(raw_pairs))
        idx_i = pairs_array[:, 0]
        idx_j = pairs_array[:, 1]
        
        vecs = coords[idx_i] - coords[idx_j]
        dists = np.linalg.norm(vecs, axis=1)
        
        types_i = np.array([elem_to_int[s] for s in symbols[idx_i]])
        types_j = np.array([elem_to_int[s] for s in symbols[idx_j]])
        
        allowed_dists = cutoff_matrix[types_i, types_j]
        
        valid_mask = (dists <= allowed_dists) & (dists > 0.1) & (allowed_dists > 0.0)
        final_pairs = pairs_array[valid_mask]
        
        existing_bonds = set()
        for b in packed_mol.bonds:
            i, j = sorted((b.atoms[0].index, b.atoms[1].index))
            existing_bonds.add((i, j))
            
        for i, j in final_pairs:
            id_i, id_j = all_atoms[i].index, all_atoms[j].index
            if id_i > id_j: id_i, id_j = id_j, id_i
            if (id_i, id_j) in existing_bonds: continue
            try:
                packed_mol.add_bond(1, all_atoms[i], all_atoms[j])
                existing_bonds.add((id_i, id_j))
            except: pass

    # --- 4. SELECT CLUSTER ---
    target_candidates = []
    for a in packed_mol.atoms:
        if a.atomic_symbol == target_element.title():
            is_match = False
            if a.label == target_atom_label: is_match = True
            elif a.label.startswith(target_atom_label + "_"): is_match = True
            elif a.label.startswith(target_atom_label + " "): is_match = True
            elif re.match(rf"^{target_atom_label}[A-Z'?]*$", a.label): is_match = True 
            if is_match: target_candidates.append(a)
    
    if not target_candidates:
        target_candidates = [a for a in packed_mol.atoms if a.atomic_symbol == target_element.title()]
    if not target_candidates: return None, None, False, 0

    center_coords_box = coords.mean(axis=0)
    start_metal = min(target_candidates, key=lambda a: np.linalg.norm(np.array(a.coordinates) - center_coords_box))
    
    # --- 5. CALCULATE NUCLEARITY ---
    nuc_visited = {start_metal.index}
    nuc_queue = [start_metal]
    found_metal_atoms = []
    if start_metal.is_metal: found_metal_atoms.append(start_metal)

    while nuc_queue:
        curr = nuc_queue.pop(0)
        for n in curr.neighbours:
            if n.index not in nuc_visited:
                nuc_visited.add(n.index)
                nuc_queue.append(n)
                if n.is_metal: found_metal_atoms.append(n)
    
    unique_sites = []
    for atom in found_metal_atoms:
        is_duplicate = False
        atom_coords = np.array([atom.coordinates.x, atom.coordinates.y, atom.coordinates.z])
        for site_coords in unique_sites:
            d = np.linalg.norm(atom_coords - site_coords)
            if d < 1.5: is_duplicate = True; break
        if not is_duplicate: unique_sites.append(atom_coords)
    true_metal_count = len(unique_sites)

    # --- 6. EXTRACT & CLEAN ---
    selected_indices = set()
    index_to_layer = {}

    if EXTRACTION_METHOD == "Geometric":
        # === IMPROVED GEOMETRIC: Iterative Dragnet ===
        
        # Queue for metals: (Atom Object, Layer Number)
        metal_queue = deque([(start_metal, 1)]) 
        visited_metal_indices = {start_metal.index}
        
        selected_indices.add(start_metal.index)
        index_to_layer[start_metal.index] = 1

        while metal_queue:
            curr_metal, curr_layer = metal_queue.popleft()
            curr_coords = np.array([curr_metal.coordinates.x, curr_metal.coordinates.y, curr_metal.coordinates.z])
            
            # --- DYNAMIC RADIUS SELECTION ---
            # 1: 3.8, 2: 3.8, 3: 7.2, 4+: 7.2
            if curr_layer <= 2:
                current_radius = Geometric_radius
            else:
                current_radius = Geometric_radius * 2

            # A. Cast the Dragnet
            nearby_indices = tree.query_ball_point(curr_coords, r=current_radius)
            
            for idx in nearby_indices:
                atom = all_atoms[idx]

                # --- EXCLUSION CHECK (Li, Na, K, etc) ---
                if atom.atomic_symbol in excluded_metals:
                    continue
                
                # B. Logic for METALS
                if atom.is_metal:
                    is_metalloligand = atom.atomic_symbol in metalloligands
                    next_layer = curr_layer if is_metalloligand else curr_layer + 1
                    
                    if next_layer <= num_metal_layers:
                        selected_indices.add(idx)
                        index_to_layer[idx] = next_layer
                        
                        if idx not in visited_metal_indices:
                            visited_metal_indices.add(idx)
                            metal_queue.append((atom, next_layer))
                            
                # C. Logic for LIGANDS
                else:
                    selected_indices.add(idx)
                    index_to_layer[idx] = curr_layer

        # === LIGAND HEALING (Topological Completion) ===
        bfs_queue = deque(list(selected_indices))
        
        while bfs_queue:
            curr_idx = bfs_queue.popleft()
            curr_atom = all_atoms[curr_idx]
            
            for n in curr_atom.neighbours:
                if n.index not in selected_indices:
                    # Excluded metals generally return True for is_metal, so they are excluded here automatically.
                    if not n.is_metal:
                        selected_indices.add(n.index)
                        index_to_layer[n.index] = index_to_layer.get(curr_idx, 1)
                        bfs_queue.append(n.index)
    else:
        # === TOPOLOGICAL ===
        queue = [(start_metal, 1)]
        visited = {start_metal.index}
        selected_indices.add(start_metal.index)
        index_to_layer[start_metal.index] = 1
        
        while queue:
            curr_atom, layer = queue.pop(0)
            for neighbor in curr_atom.neighbours:
                if neighbor.index in visited: continue
                
                # --- EXCLUSION CHECK ---
                if neighbor.atomic_symbol in excluded_metals: continue
                
                # --- GUARD: DOUBLE CHECK ---
                # Even with strict matrix, if a bond existed in CSD that slipped through
                # because it was < 5.0 but > 2.8, we catch it here.
                # Calculate distance dynamically
                d_real = distance(curr_atom, neighbor)
                
                # If both are Non-Metals and distance > 2.8, STOP.
                if not curr_atom.is_metal and not neighbor.is_metal:
                    if d_real > LIMIT_NM_NM: continue
                
                # If one is Metal and dist > 3.5, STOP
                elif (curr_atom.is_metal != neighbor.is_metal):
                    if d_real > LIMIT_M_NM: continue

                # If H...NonMetal, STOP (The H-bond guard)
                if curr_atom.atomic_symbol == 'H' and not neighbor.is_metal: continue

                is_metal = neighbor.is_metal
                sym = neighbor.atomic_symbol
                next_layer = layer
                if is_metal:
                    if sym in active_metalloligands: next_layer = layer
                    else: next_layer = layer + 1
                
                if next_layer <= num_metal_layers:
                    visited.add(neighbor.index)
                    selected_indices.add(neighbor.index)
                    index_to_layer[neighbor.index] = next_layer
                    
                    if is_metal:
                        if sym in active_metalloligands: queue.append((neighbor, next_layer))
                        elif next_layer < num_metal_layers: queue.append((neighbor, next_layer))
                    else:
                        queue.append((neighbor, next_layer))

    # --- COMMON CLEANUP ---
    selected_list = list(selected_indices)
    selected_coords = coords[selected_list]
    overlaps = tree.query_ball_point(selected_coords, r=1.0)
    
    final_indices = set(selected_indices)
    for i, nearby_indices in enumerate(overlaps):
        original_idx = selected_list[i]
        base_layer = index_to_layer.get(original_idx, 1)
        for nearby_idx in nearby_indices:
            if all_atoms[nearby_idx].atomic_symbol in excluded_metals: continue
            if nearby_idx not in final_indices:
                final_indices.add(nearby_idx)
                index_to_layer[nearby_idx] = base_layer

    clean_mol = molecule.Molecule(identifier=f"{entry.identifier}_{target_atom_label}")
    atom_to_layer_map = {}
    old_idx_to_new_atom = {}
    
    for idx in final_indices:
        original = all_atoms[idx]
        new_atom = clean_mol.add_atom(original)
        old_idx_to_new_atom[idx] = new_atom
        layer = index_to_layer.get(idx, 1)
        atom_to_layer_map[new_atom] = layer
        
    # We check the ORIGINAL connectivity in 'all_atoms' (the packed cell)
    # If an atom was connected to >1 metal in the original crystal, it IS bridging,
    # even if we stripped the other metals in the extraction.
    proven_bridging_indices = set()
    for idx in final_indices:
        original_atom = all_atoms[idx]
        if original_atom.is_metal: continue
        
        # Count metal neighbors in the full PARENT structure
        parent_metal_neighbors = 0
        for n in original_atom.neighbours:
            if n.is_metal: 
                parent_metal_neighbors += 1
        
        if parent_metal_neighbors > 1:
            # It is bridging! Record the index in the NEW clean_mol
            new_atom = old_idx_to_new_atom[idx]
            proven_bridging_indices.add(new_atom.index)
            
    subset_indices = list(final_indices)
    subset_set = set(subset_indices)
    
    for old_idx in subset_indices:
        original = all_atoms[old_idx]
        new_a1 = old_idx_to_new_atom[old_idx]
        for bond in original.bonds:
            other = bond.atoms[0] if bond.atoms[1].index == original.index else bond.atoms[1]
            if other.index in subset_set and other.index > old_idx:
                new_a2 = old_idx_to_new_atom[other.index]
                try: clean_mol.add_bond(bond.bond_type, new_a1, new_a2)
                except RuntimeError: clean_mol.add_bond(1, new_a1, new_a2)

    # --- FIXED MANUAL COMMAND BLOCK ---
    if manual_commands:
        # 1. Create a stable mapping of the current indices to atom objects
        id_map = {a.index: a for a in clean_mol.atoms}
        
        all_atoms_to_remove = []

        # A. Collect Explicit Deletions
        indices_to_del = manual_commands.get('delete', [])
        for idx in indices_to_del:
            if idx in id_map:
                all_atoms_to_remove.append(id_map[idx])

        # B. Handle Averaging
        avg_groups = manual_commands.get('average', [])
        for group in avg_groups:
            target_atoms = [id_map.get(idx) for idx in group if idx in id_map]
            
            if len(target_atoms) > 1:
                leader = target_atoms[0]
                followers = target_atoms[1:]

                # Calculate weighted coordinates
                total_occ = sum(float(getattr(a, 'occupancy', 1.0) or 1.0) for a in target_atoms)
                
                if total_occ > 0:
                    avg_x = sum(a.coordinates.x * float(getattr(a, 'occupancy', 1.0) or 1.0) for a in target_atoms) / total_occ
                    avg_y = sum(a.coordinates.y * float(getattr(a, 'occupancy', 1.0) or 1.0) for a in target_atoms) / total_occ
                    avg_z = sum(a.coordinates.z * float(getattr(a, 'occupancy', 1.0) or 1.0) for a in target_atoms) / total_occ
                    
                    # Update leader's position (this is allowed)
                    leader.coordinates = molecule.Coordinates(avg_x, avg_y, avg_z)
                    
                    # Mark all other atoms in the group for removal
                    all_atoms_to_remove.extend(followers)

        # C. Single-Pass Removal to prevent index shifting
        if all_atoms_to_remove:
            clean_mol.remove_atoms(list(set(all_atoms_to_remove)))
        
        found_disorder = True
    else:
        clean_mol, found_disorder = clean_small_cluster(clean_mol, atom_to_layer_map, 
                                                        perform_correction=perform_correction,
                                                        disorder_method=disorder_resolve_method)
        

    # --- FINAL VALIDATION ---
    new_center = old_idx_to_new_atom.get(start_metal.index)
    
    if not new_center: 
        return None, None, found_disorder, 0, set(), True # Failed to find center
    
    validation_failed = False
    # We run validation but DO NOT return None if it fails.
    # Instead, we just mark the flag.
    is_valid = validate_geometry(clean_mol, new_center, check_rdkit=True, emit=emit)
    if not is_valid: 
         validation_failed = True
        
    return clean_mol, new_center, found_disorder, true_metal_count, proven_bridging_indices, validation_failed


def run_pipeline(
    config: Config,
    emit: Callable[[str], None] = default_emit,
    display_views: bool = True,
    on_visualization: Optional[Callable[[str, int], None]] = None,
) -> Result:
    all_results_list: List[Dict[str, object]] = []
    visualizations_html: List[str] = []
    summary_reports: List[str] = []
    artifacts: Dict[str, Optional[object]] = {
        "run_summary_txt": "run_summary.txt",
        "pkl_paths": [],
        "csv_paths": [],
    }
    abstract_cache: Dict[str, Optional[str]] = {}

    def _iter_spool(spool_path: str):
        with open(spool_path, "rb") as spool_handle:
            while True:
                try:
                    yield pickle.load(spool_handle)
                except EOFError:
                    break
                except pickle.UnpicklingError:
                    break

    def _read_spool(spool_path: str) -> List[Dict[str, object]]:
        return list(_iter_spool(spool_path))

    def _filter_timeout_rows(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
        filtered = []
        for row in rows:
            timeout_site = _parse_timeout_site(row.get("CSD ID"))
            if timeout_site:
                continue
            filtered.append(row)
        return filtered

    def _has_unresolved_timeouts(rows: List[Dict[str, object]]) -> bool:
        finished_sites = set()
        timeout_sites = set()
        for row in rows:
            timeout_site = _parse_timeout_site(row.get("CSD ID"))
            if timeout_site:
                timeout_sites.add(timeout_site)
                continue
            csd_id = row.get("CSD ID")
            site_label = row.get("Site Label")
            if csd_id and site_label:
                finished_sites.add(f"{csd_id}_{site_label}")
        return any(site not in finished_sites for site in timeout_sites)

    def _load_resume_state(spool_path: str, retry_timeouts: bool):
        resume_seen = set()
        resume_stats = {"Text Mining (Metadata)": 0, "Gemini LLM (Abstract)": 0, "BVS (Geometry)": 0}
        resume_processed = 0
        try:
            for row in _iter_spool(spool_path):
                timeout_site = _parse_timeout_site(row.get("CSD ID"))
                if timeout_site:
                    if not retry_timeouts and timeout_site != "timeout":
                        resume_processed += 1
                        resume_seen.add(timeout_site)
                    continue
                resume_processed += 1
                csd_id = row.get("CSD ID")
                site_label = row.get("Site Label")
                if csd_id and site_label:
                    resume_seen.add(f"{csd_id}_{site_label}")
                method_used = str(row.get("OS_Method", ""))
                if "Text" in method_used:
                    resume_stats["Text Mining (Metadata)"] += 1
                elif "Gemini" in method_used:
                    resume_stats["Gemini LLM (Abstract)"] += 1
                elif "BVS" in method_used:
                    resume_stats["BVS (Geometry)"] += 1
        except Exception as exc:
            emit(f"Warning: could not read resume spool {spool_path}: {exc}")
        return resume_seen, resume_processed, resume_stats

    for TARGET_METAL in config.TARGET_METAL_list:
        hits = []
        suffix = "specific_list" if config.target_csd_ids else "search"
        pkl_name = f"{TARGET_METAL}_layers{config.num_metal_layers}_{suffix}_results.pkl"
        csv_name = f"{TARGET_METAL}_layers{config.num_metal_layers}_{suffix}_results.csv"
        spool_name = f".{TARGET_METAL}_layers{config.num_metal_layers}_{suffix}_results.spool.pkl"
        spool_handle = None
        spool_count = 0

        # --- BRANCH 1: LIST OF CSD IDS PROVIDED ---
        if config.target_csd_ids and len(config.target_csd_ids) > 0:
            emit(f"Initializing Analysis for specific list of {len(config.target_csd_ids)} entries...")
            entry_reader = io.EntryReader("CSD")
            found_count = 0
            for refcode in config.target_csd_ids:
                try:
                    clean_ref = str(refcode).strip()
                    entry = entry_reader.entry(clean_ref)
                    if entry:
                        hits.append(SimpleHit(entry))
                        found_count += 1
                    else:
                        emit(f"Warning: Refcode {clean_ref} not found in CSD.")
                except Exception as e:
                    emit(f"Error retrieving {refcode}: {e}")
            emit(f"Successfully retrieved {found_count} entries.")

        # --- BRANCH 2: SEARCH BY ELEMENT ---
        else:
            emit(f"Initializing CSD Search for {TARGET_METAL}...")
            metal_search = search.SubstructureSearch()
            metal_search.add_substructure(search.SMARTSSubstructure(f"[{TARGET_METAL}]"))
            hits = metal_search.search()
            emit(f"Found {len(hits)} total structures containing {TARGET_METAL}.")

        if not hits:
            emit("No structures found to process.")
            continue

        # --- SETUP STATISTICS ---
        stats_counts = {"Text Mining (Metadata)": 0, "Gemini LLM (Abstract)": 0, "BVS (Geometry)": 0}
        skipped_r_factor = 0
        skipped_polymeric = 0
        skipped_powder = 0
        skipped_all_disorder = 0
        skipped_primary_disorder = 0
        skipped_extraction_fail = 0
        skipped_os_mismatch = 0
        skipped_timeout = 0
        corrected_disorder_count = 0
        correction_failed_count = 0

        processed_count = 0
        visualized_count = 0

        seen_site_ids = set()
        if os.path.exists(spool_name):
            resume_seen, resume_processed, resume_stats = _load_resume_state(
                spool_name,
                config.RETRY_TIMEOUTS,
            )
            if resume_processed > 0:
                seen_site_ids = resume_seen
                processed_count = resume_processed
                for key in stats_counts:
                    stats_counts[key] += resume_stats.get(key, 0)
                emit(f"Resuming from {spool_name}: {resume_processed} results already collected.")

        emit("\n--- PROCESSING STRUCTURES ---\n")

        for hit in hits:
            if processed_count >= config.PROCESS_LIMIT:
                break

            # --- METADATA FILTERS ---
            if hit.entry.r_factor is None or hit.entry.r_factor >= config.R_FACTOR_LIMIT:
                skipped_r_factor += 1
                continue

            if config.FILTER_POLYMERIC and hit.entry.is_polymeric:
                skipped_polymeric += 1
                continue

            if config.FILTER_POWDER:
                is_powder = False
                if hit.entry.synonyms:
                    for s in hit.entry.synonyms:
                        if "powder" in s.lower():
                            is_powder = True
                            break
                if is_powder:
                    skipped_powder += 1
                    continue

            if config.FILTER_ALL_DISORDER and hit.entry.has_disorder:
                skipped_all_disorder += 1
                continue

            # --- MULTI-SITE ITERATION LOGIC ---
            try:
                entry_atoms = hit.entry.crystal.asymmetric_unit_molecule.atoms
            except AttributeError:
                entry_atoms = hit.entry.molecule.atoms

            target_atoms_in_entry = [a for a in entry_atoms if a.atomic_symbol == TARGET_METAL.title()]

            if not target_atoms_in_entry:
                continue

            unique_targets = []
            seen_base_labels_in_entry = set()

            for atom in target_atoms_in_entry:
                match = re.match(r"([A-Za-z]{1,2}\d+)", atom.label)
                if match:
                    base_label = match.group(1)
                else:
                    base_label = atom.label

                if base_label not in seen_base_labels_in_entry:
                    seen_base_labels_in_entry.add(base_label)
                    unique_targets.append((atom, base_label))

            for target_atom, clean_base_label in unique_targets:
                unique_site_id = f"{hit.identifier}_{clean_base_label}"
                if unique_site_id in seen_site_ids:
                    continue
                seen_site_ids.add(unique_site_id)

                use_timeout = config.SITE_TIMEOUT_SECONDS and config.SITE_TIMEOUT_SECONDS > 0
                if use_timeout:
                    payload = {
                        "csd_id": hit.identifier,
                        "target_metal": TARGET_METAL,
                        "clean_base_label": clean_base_label,
                        "config": config,
                    }
                    timeout_result = _run_site_with_timeout(payload, config.SITE_TIMEOUT_SECONDS)
                    status = timeout_result.get("status")
                    if status == "timeout":
                        skipped_timeout += 1
                        emit(
                            f"Timeout: {unique_site_id} exceeded {config.SITE_TIMEOUT_SECONDS} seconds. Skipping."
                        )
                        timeout_value = _build_timeout_value(hit.identifier, clean_base_label)
                        timeout_row = {key: timeout_value for key in RESULT_COLUMNS}
                        if spool_handle is None:
                            spool_handle = open(spool_name, "ab")
                        pickle.dump(timeout_row, spool_handle, protocol=pickle.HIGHEST_PROTOCOL)
                        spool_count += 1
                        if spool_count % 50 == 0:
                            spool_handle.flush()
                        continue
                    if status == "fail_correction":
                        correction_failed_count += 1
                        continue
                    if status == "fail_extraction":
                        skipped_extraction_fail += 1
                        continue
                    if status == "skip_primary_disorder":
                        skipped_primary_disorder += 1
                        continue
                    if status == "skip_os_mismatch":
                        skipped_os_mismatch += 1
                        continue
                    if status == "entry_missing":
                        emit(f"Processing failed for {unique_site_id}: entry not found in CSD.")
                        skipped_extraction_fail += 1
                        continue
                    if status == "error":
                        emit(f"Processing failed for {unique_site_id}: {timeout_result.get('error')}")
                        skipped_extraction_fail += 1
                        continue

                    result_row = timeout_result.get("result_row")
                    if not result_row:
                        emit(f"Processing failed for {unique_site_id}: no result produced.")
                        skipped_extraction_fail += 1
                        continue

                    if timeout_result.get("has_primary_disorder") and config.CORRECT_PRIMARY_DISORDER:
                        corrected_disorder_count += 1

                    processed_count += 1
                    validation_failed = result_row.get("VALIDATION_FAILED", False)
                    emit(
                        f"{processed_count}. {unique_site_id} [OS: {result_row.get('Oxidation State')}] "
                        f"[Valid: {not validation_failed}]"
                    )

                    method_used = timeout_result.get("method_used", "")
                    if "Text" in method_used:
                        stats_counts["Text Mining (Metadata)"] += 1
                    elif "Gemini" in method_used:
                        stats_counts["Gemini LLM (Abstract)"] += 1
                    elif "BVS" in method_used:
                        stats_counts["BVS (Geometry)"] += 1

                    if spool_handle is None:
                        spool_handle = open(spool_name, "ab")
                    pickle.dump(result_row, spool_handle, protocol=pickle.HIGHEST_PROTOCOL)
                    spool_count += 1
                    if spool_count % 50 == 0:
                        spool_handle.flush()

                    if config.VISUALIZE and visualized_count < config.VISUALIZATION_LIMIT:
                        AN_n_xyz_block = result_row.get("AN_n_XYZ Coordinates", "")
                        view = py3Dmol.view(width=400, height=400)
                        view.addModel(AN_n_xyz_block, "xyz")
                        view.setStyle({"stick": {}, "sphere": {"scale": 0.2}})
                        view.addStyle({"elem": TARGET_METAL}, {"sphere": {"scale": 0.5, "color": "#00FF00"}})
                        view.zoomTo()
                        html = view._make_html()
                        visualizations_html.append(html)
                        if on_visualization:
                            on_visualization(html, visualized_count + 1)
                        if display_views:
                            view.show()
                        visualized_count += 1
                    continue

                # === EXPANSION, CLEANING, CORRECTION ===
                try:
                    core_result = process_structure_core(
                        hit.entry,
                        TARGET_METAL,
                        target_atom_label=clean_base_label,
                        num_metal_layers=config.num_metal_layers,
                        perform_correction=config.CORRECT_PRIMARY_DISORDER,
                        EXTRACTION_METHOD=config.EXTRACTION_METHOD,
                        disorder_resolve_method=config.DISORDER_RESOLVE_METHOD,
                        Extraction_Cut_off_distances=config.Extraction_Cut_off_distances,
                        Geometric_radius=config.Geometric_radius,
                        metalloligands=config.metalloligands,
                        emit=emit,
                    )
                    mol, center_metal, has_primary_disorder, nuclearity_count, bridging_indices, validation_failed = core_result

                    if mol is None:
                        if has_primary_disorder and config.CORRECT_PRIMARY_DISORDER:
                            correction_failed_count += 1
                        else:
                            skipped_extraction_fail += 1
                        continue

                    if has_primary_disorder:
                        if config.FILTER_PRIMARY_DISORDER:
                            skipped_primary_disorder += 1
                            continue
                        if config.CORRECT_PRIMARY_DISORDER:
                            corrected_disorder_count += 1

                except Exception as e:
                    emit(f"Processing failed for {unique_site_id}: {e}")
                    skipped_extraction_fail += 1
                    continue

                # === ANALYSIS ===
                name_str = hit.entry.chemical_name if hit.entry.chemical_name else ""
                form_str = hit.entry.formula if hit.entry.formula else ""
                combined_text = f"{name_str} {form_str}"

                # 1. Ligand Analysis
                AN_n_ligand_data, has_bridging, o_map_n, virtual_hs = analyze_ligands(
                    mol,
                    center_metal,
                    combined_text,
                    bridging_indices,
                    Hydrogen_Addition_method=config.Hydrogen_Addition_method,
                )

                # --- POST-HYDROGEN VALIDATION ---
                is_valid = validate_geometry(
                    mol,
                    center_metal,
                    virtual_hydrogens=virtual_hs,
                    check_rdkit=True,
                    emit=emit,
                )
                validation_failed = not is_valid

                # 2. XYZ Generation
                AN_n_xyz_block = generate_enriched_xyz_string(mol, virtual_hs)

                # 3. Oxidation States Calculation
                bvs_val = calculate_bvs(center_metal, center_metal.neighbours)
                bvs_int = round(bvs_val)
                final_state = "?"
                source = "Unknown"
                source_content = ""
                method_used = ""

                candidate_states_list = []

                # --- ABSTRACT FETCHING LOGIC ---
                abstract_text = "N/A"

                should_fetch = config.GET_ABSTRACT or HAS_GEMINI

                if should_fetch and hit.entry.publication and hit.entry.publication.doi:
                    doi = str(hit.entry.publication.doi).strip()
                    doi_key = doi.lower()
                    if doi_key in abstract_cache:
                        cached_abstract = abstract_cache[doi_key]
                    else:
                        cached_abstract = get_abstract_from_openalex(doi_key)
                        abstract_cache[doi_key] = cached_abstract
                    if cached_abstract:
                        abstract_text = cached_abstract

                # --- MINING ---
                text_states, combined_text_os = extract_text_oxidation_state(hit.entry, TARGET_METAL)
                if text_states:
                    candidate_states_list.extend(text_states)
                    source = "Text Mining"
                    method_used = "Text Mining"

                # Only run Gemini if we successfully got an abstract
                elif HAS_GEMINI and abstract_text != "N/A":
                    llm_states = extract_oxidation_state_with_gemini(abstract_text, TARGET_METAL)
                    if llm_states:
                        candidate_states_list.extend(llm_states)
                        source = "Gemini LLM"
                        method_used = "Gemini LLM"

                if len(candidate_states_list) > 1:
                    final_state = resolve_mixed_valence(center_metal, candidate_states_list, mol)
                    method_used += " (Resolved by Geometry)"
                elif len(candidate_states_list) == 1:
                    final_state = f"+{candidate_states_list[0]}"
                    source = source or "Text Mining"
                    method_used = method_used or "Text Mining"
                else:
                    final_state = f"+{bvs_int}"
                    source = "BVS"
                    method_used = "BVS (Geometry)"

                # --- FILTER BY USER PROVIDED STATES ---
                should_keep = False
                if not config.OXIDATION_STATES_FILTER or config.OXIDATION_STATES_FILTER == ["all"]:
                    should_keep = True
                else:
                    for allowed_state in config.OXIDATION_STATES_FILTER:
                        if allowed_state.lower() == "actinyl":
                            if "actinyl" in (source_content or combined_text).lower():
                                should_keep = True
                                break
                        else:
                            try:
                                allowed_int = int(str(allowed_state).replace("+", ""))
                                if final_state.lstrip("+") == str(allowed_int):
                                    should_keep = True
                                    break
                            except ValueError:
                                continue

                if not should_keep:
                    skipped_os_mismatch += 1
                    continue

                processed_count += 1
                emit(f"{processed_count}. {unique_site_id} [OS: {final_state}] [Valid: {not validation_failed}]")

                if "Text" in method_used:
                    stats_counts["Text Mining (Metadata)"] += 1
                elif "Gemini" in method_used:
                    stats_counts["Gemini LLM (Abstract)"] += 1
                elif "BVS" in method_used:
                    stats_counts["BVS (Geometry)"] += 1

                structure_type = get_nuclearity_label(nuclearity_count)
                coordination_number = sum(len(x["Connecting_atom"]) for x in AN_n_ligand_data)
                total_ligand_charge = sum(item["Charge"] * item["number"] for item in AN_n_ligand_data)
                lig_count = sum(x["number"] for x in AN_n_ligand_data)
                has_rdkit_issue = any(item.get("HAD_RDKIT_ISSUE", False) for item in AN_n_ligand_data)

                doi_str = hit.entry.publication.doi if hit.entry.publication else "N/A"

                # --- EXTRACT UNIT CELL INFO ---
                uc_dict = {}
                if hit.entry.crystal:
                    uc_dict = {
                        "a": hit.entry.crystal.cell_lengths.a,
                        "b": hit.entry.crystal.cell_lengths.b,
                        "c": hit.entry.crystal.cell_lengths.c,
                        "alpha": hit.entry.crystal.cell_angles.alpha,
                        "beta": hit.entry.crystal.cell_angles.beta,
                        "gamma": hit.entry.crystal.cell_angles.gamma,
                    }

                result_row = {
                    "CSD ID": hit.identifier,
                    "Site Label": clean_base_label,
                    "Metal": TARGET_METAL,
                    "Oxidation State": final_state,
                    "Chemical Name": name_str,
                    "Unit_cell": uc_dict,
                    "Abstract": abstract_text,
                    "Structure Type": structure_type,
                    "Coordination Number": coordination_number,
                    "Ligand_count": lig_count,
                    "Ligand_Info": AN_n_ligand_data,
                    "VALIDATION_FAILED": validation_failed,
                    "HAD_RDKIT_ISSUE": has_rdkit_issue,
                    "Bridging_Ligand": has_bridging,
                    "Total_Ligand_Charge": total_ligand_charge,
                    "OS_Method": method_used,
                    "Hydrogen Method": config.Hydrogen_Addition_method,
                    "OS_source": source,
                    "XYZ Coordinates": AN_n_xyz_block,
                    "Had_prismary_sphere_disorder": has_primary_disorder,
                    "DOI": doi_str,
                }
                if spool_handle is None:
                    spool_handle = open(spool_name, "ab")
                pickle.dump(result_row, spool_handle, protocol=pickle.HIGHEST_PROTOCOL)
                spool_count += 1
                if spool_count % 50 == 0:
                    spool_handle.flush()

                if config.VISUALIZE and visualized_count < config.VISUALIZATION_LIMIT:
                    view = py3Dmol.view(width=400, height=400)
                    view.addModel(AN_n_xyz_block, "xyz")
                    view.setStyle({"stick": {}, "sphere": {"scale": 0.2}})
                    view.addStyle({"elem": TARGET_METAL}, {"sphere": {"scale": 0.5, "color": "#00FF00"}})
                    view.zoomTo()
                    html = view._make_html()
                    visualizations_html.append(html)
                    if on_visualization:
                        on_visualization(html, visualized_count + 1)
                    if display_views:
                        view.show()
                    visualized_count += 1

        # --- FINAL SUMMARY ---
        summary_report = (
            "\n"
            + "=" * 50
            + "\n"
            + f" SUMMARY STATISTICS ({TARGET_METAL} Processed {processed_count} unique sites)\n"
            + "=" * 50
            + "\n"
            + f" Filter Mode (Oxidation States): {config.OXIDATION_STATES_FILTER}\n"
            + f" Metal Layers (Depth): {config.num_metal_layers}\n"
            + f" Total Hits: {len(hits)}\n"
            + "-" * 50
            + "\n"
            + f" Skipped (R-Factor >= {config.R_FACTOR_LIMIT}% or None): {skipped_r_factor}\n"
            + f" Skipped (Polymeric): {skipped_polymeric}\n"
            + f" Skipped (Powder Data): {skipped_powder}\n"
            + f" Skipped (Has Global Disorder): {skipped_all_disorder}\n"
            + f" Skipped (Has Primary Disorder): {skipped_primary_disorder}\n"
            + f" Skipped (Extraction/Validation Failed): {skipped_extraction_fail}\n"
            + f" Skipped (Oxidation State Filter Mismatch): {skipped_os_mismatch}\n"
            + f" Skipped (Correction Failed): {correction_failed_count}\n"
            + f" Skipped (Timeout): {skipped_timeout}\n"
            + f" Structures Corrected Successfully: {corrected_disorder_count}\n"
            + "-" * 50
            + "\n"
        )
        emit(summary_report)
        summary_reports.append(summary_report)
        with open("run_summary.txt", "a") as f:
            import datetime

            timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            f.write(f"Run Timestamp: {timestamp}\n{summary_report}\n\n")

        if spool_handle is not None:
            spool_handle.close()
        raw_spool_rows = _read_spool(spool_name) if os.path.exists(spool_name) else []
        results_list = _filter_timeout_rows(raw_spool_rows) if raw_spool_rows else []
        if results_list:
            try:
                df = pd.DataFrame(results_list)
                df.to_pickle(pkl_name)
                df.to_csv(csv_name, index=False)
                artifacts["pkl_paths"].append(pkl_name)
                artifacts["csv_paths"].append(csv_name)
                emit(f" ✅ Data exported to {pkl_name} and {csv_name}")
            except Exception as e:
                emit(f" ❌ Export failed: {e}")
        else:
            emit(" No structures processed.")
        if raw_spool_rows and not _has_unresolved_timeouts(raw_spool_rows):
            if os.path.exists(spool_name):
                os.remove(spool_name)
        emit("=" * 10)

        if config.Edit_manual and processed_count > 0:
            emit("Running Manual Refinment")
            from .Manual import manual_refinement_step
            manual_refinement_step(
                TARGET_METAL,
                config.num_metal_layers,
                hits,
                config.Hydrogen_Addition_method,
                config.Extraction_Cut_off_distances,
                pkl_name=pkl_name,
                csv_name=csv_name,
                EXTRACTION_METHOD=config.EXTRACTION_METHOD,
                Geometric_radius=config.Geometric_radius,
            )

    df = None
    if artifacts["pkl_paths"]:
        frames = []
        for pkl_path in artifacts["pkl_paths"]:
            try:
                frames.append(pd.read_pickle(pkl_path))
            except Exception as e:
                emit(f" ❌ Aggregate read failed for {pkl_path}: {e}")
        if frames:
            df = pd.concat(frames, ignore_index=True)
            all_results_list = df.to_dict(orient="records")
    return Result(
        results=all_results_list,
        dataframe=df,
        summary_report="".join(summary_reports),
        artifacts=artifacts,
        visualizations_html=visualizations_html,
    )
