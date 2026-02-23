import math
import numpy as np

from .constants import ELEMENT_NAMES, actinyl_thresholds
from .utils import normalize, get_vec, distance
try:
    from rdkit import Chem
    from rdkit import Geometry
    from rdkit.Chem import AllChem, rdDetermineBonds
    from rdkit import RDLogger
    from rdkit import rdBase
    RDLogger.EnableLog('rdApp.error')
    rdBase.LogToPythonStderr()
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

def calculate_virtual_hydrogens(atom, neighbor_atoms, geometry_type="tetrahedral", target_H_count=1):
    """
    Calculates coordinates for missing Hydrogens based on geometry.
    Supports sp (linear), sp2 (planar), sp3 (tetrahedral).
    """
    atom_coords = np.array([atom.coordinates.x, atom.coordinates.y, atom.coordinates.z])
    new_h_coords = []
    
    # Standard Bond Lengths
    H_DIST_SP3 = 1.09
    H_DIST_SP2 = 1.08
    H_DIST_SP  = 1.06
    
    # Vectors to ALL existing neighbors (Heavy + existing H)
    vectors = [normalize(get_vec(atom, n)) for n in neighbor_atoms]
    n_neigh = len(vectors)

    # --- 1. SP LINEAR (Alkynes, Nitriles) ---
    if geometry_type == "linear" or geometry_type == "sp_carbon":
        if n_neigh == 1:
            # Point 180 deg away
            u = -vectors[0]
            new_h_coords.append(tuple(atom_coords + u * H_DIST_SP))

    # --- 2. SP2 PLANAR (Aromatic / Alkenes) ---
    elif geometry_type == "sp2_carbon" or geometry_type == "aromatic":
        if n_neigh == 2:
            # Ideal trigonal planar missing 1 atom
            v1, v2 = vectors[0], vectors[1]
            sum_v = v1 + v2
            if np.linalg.norm(sum_v) < 0.1: # Error case: Linear neighbors
                u = normalize(np.cross(v1, np.array([0,0,1])))
            else:
                u = normalize(-sum_v) # Bisector
            new_h_coords.append(tuple(atom_coords + u * H_DIST_SP2))
            
        elif n_neigh == 1:
            # Terminal double bond (needs 2 H)
            v1 = vectors[0]
            arb = np.array([0,0,1]) if abs(v1[2]) < 0.9 else np.array([0,1,0])
            plane_norm = normalize(np.cross(v1, arb))
            u = normalize(np.cross(plane_norm, v1))
            
            # Add at +/- 120 deg
            h1 = atom_coords + (-0.5 * v1 + 0.866 * u) * H_DIST_SP2
            h2 = atom_coords + (-0.5 * v1 - 0.866 * u) * H_DIST_SP2
            new_h_coords = [tuple(h1), tuple(h2)]

    # --- 3. SP3 TETRAHEDRAL (Aliphatic) ---
    elif geometry_type == "sp3_carbon" or geometry_type == "tetrahedral":
        if n_neigh == 3:
            # Tertiary (Tripod missing top)
            sum_v = vectors[0] + vectors[1] + vectors[2]
            u = normalize(-sum_v)
            new_h_coords.append(tuple(atom_coords + u * H_DIST_SP3))
            
        elif n_neigh == 2:
            # Secondary (V-shape missing 2 legs)
            v1, v2 = vectors[0], vectors[1]
            bisector = normalize(v1 + v2)
            perp = normalize(np.cross(v1, v2))
            u_main = -bisector 
            # Tetrahedral projection
            h1 = atom_coords + (u_main * 0.5 + perp * 0.866) * H_DIST_SP3
            h2 = atom_coords + (u_main * 0.5 - perp * 0.866) * H_DIST_SP3
            new_h_coords = [tuple(h1), tuple(h2)]
            
        elif n_neigh == 1:
            # Primary (Line missing 3 legs)
            v1 = vectors[0]
            u = -v1 
            arb = np.array([0,0,1]) if abs(u[2]) < 0.9 else np.array([0,1,0])
            v_perp1 = normalize(np.cross(u, arb))
            v_perp2 = normalize(np.cross(u, v_perp1))
            
            # Cone at 109.5 deg
            cos_t = 0.33333; sin_t = 0.94281
            h1 = atom_coords + (u * cos_t + v_perp1 * sin_t) * H_DIST_SP3
            v_p2 = -0.5 * v_perp1 + 0.866 * v_perp2
            v_p3 = -0.5 * v_perp1 - 0.866 * v_perp2
            h2 = atom_coords + (u * cos_t + v_p2 * sin_t) * H_DIST_SP3
            h3 = atom_coords + (u * cos_t + v_p3 * sin_t) * H_DIST_SP3
            new_h_coords = [tuple(h1), tuple(h2), tuple(h3)]

    # --- 4. WATER / HYDROXYL ---
    elif geometry_type == "water_like":
        H_DIST = 0.96
        if n_neigh == 1:
            u = -vectors[0]
            arb = np.array([0,0,1]) if abs(u[2]) < 0.9 else np.array([0,1,0])
            v1 = normalize(np.cross(u, arb))
            v2 = np.cross(u, v1)
            sin_t = math.sin(math.radians(55)); cos_t = math.cos(math.radians(55))
            h1 = atom_coords + (u * cos_t + v1 * sin_t) * H_DIST
            h2 = atom_coords + (u * cos_t - v1 * sin_t * 0.5 + v2 * sin_t * 0.866) * H_DIST
            new_h_coords = [tuple(h1), tuple(h2)]

    elif geometry_type == "hydroxyl_like":
        H_DIST = 0.96
        if n_neigh == 1:
            u = -vectors[0]
            arb = np.array([0,0,1]) if abs(u[2]) < 0.9 else np.array([0,1,0])
            perp = normalize(np.cross(u, arb))
            sin_t = math.sin(math.radians(70)); cos_t = math.cos(math.radians(70))
            h1 = atom_coords + (u * cos_t + perp * sin_t) * H_DIST
            new_h_coords = [tuple(h1)]
            
    return new_h_coords[:target_H_count]

class LigandProcessor:
    def __init__(self, mol, metal_atom, text_source="", known_bridging_indices=None, 
                 backbone_tolerance=0.3, Hydrogen_Addition_method='Geometric'):
        """
        Initializes the processor with a specific hydrogen addition method.
        Hydrogen_Addition_method options: 'None', 'Geometric', 'RDkit'
        """
        self.mol = mol
        self.metal_atom = metal_atom
        self.text_source = text_source.lower()
        self.known_bridging_indices = known_bridging_indices if known_bridging_indices else set()
        self.backbone_tolerance = backbone_tolerance 
        self.addition_method = Hydrogen_Addition_method # New parameter for the 3 choices
        self.virtual_hydrogens = [] 
        self.ligand_data = []
        self.has_bridging = False
        self.oxygen_map = {} 
        
        self._build_atom_map()
        self._process()
        
    def _build_atom_map(self):
        self.mol_atom_map = []
        for atom in self.mol.atoms:
            if atom.coordinates:
                c = np.array([atom.coordinates.x, atom.coordinates.y, atom.coordinates.z])
                self.mol_atom_map.append((atom, c))

    def _get_original_atom(self, component_atom):
        if not component_atom.coordinates: return None
        target = np.array([component_atom.coordinates.x, component_atom.coordinates.y, component_atom.coordinates.z])
        best_atom = None
        min_dist = 0.1 
        
        for atom, coords in self.mol_atom_map:
            if atom.atomic_symbol != component_atom.atomic_symbol: continue
            dist = np.linalg.norm(coords - target)
            if dist < min_dist:
                min_dist = dist
                best_atom = atom
        return best_atom

    def _is_borate_network(self, component):
        symbols = {a.atomic_symbol for a in component.atoms}
        return "B" in symbols and symbols.issubset({"B", "O", "H"})

    def _collect_borate_fragments(self, component):
        atom_by_index = {a.index: a for a in component.atoms}
        component_indices = set(atom_by_index.keys())
        fragments = []
        seen = set()
        for atom in component.atoms:
            if atom.atomic_symbol != "B":
                continue
            keep_indices = {atom.index}
            for neighbor in atom.neighbours:
                if neighbor.index not in component_indices:
                    continue
                keep_indices.add(neighbor.index)
                for h_neighbor in neighbor.neighbours:
                    if h_neighbor.atomic_symbol == "H" and h_neighbor.index in component_indices:
                        keep_indices.add(h_neighbor.index)
            key = frozenset(keep_indices)
            if key not in seen:
                seen.add(key)
                fragments.append(keep_indices)
        return fragments

    def _process_borate_component(self, component):
        for keep_indices in self._collect_borate_fragments(component):
            fragment = component.copy()
            atoms_to_remove = [a for a in fragment.atoms if a.index not in keep_indices]
            if atoms_to_remove:
                fragment.remove_atoms(atoms_to_remove)
            connecting_atoms = []
            for comp_atom in fragment.atoms:
                orig_atom = self._get_original_atom(comp_atom)
                if not orig_atom:
                    continue
                for neighbor in orig_atom.neighbours:
                    if neighbor.index == self.metal_atom.index:
                        connecting_atoms.append(orig_atom)
                        break
            if connecting_atoms:
                self._analyze_single_ligand(fragment, connecting_atoms)

    # --- UNCHANGED: _determine_geometry_and_valency ---
    def _determine_geometry_and_valency(self, atom, neighbors):
        """
        Determines target geometry and valency using Bond Valence logic for N,
        and geometric logic for C.
        """
        sym = atom.atomic_symbol.strip()
        n_neigh = len(neighbors)
        atom_c = np.array([atom.coordinates.x, atom.coordinates.y, atom.coordinates.z])

        # === HELPER: Nitrogen Bond Order Calculator ===
        def get_n_bond_order(dist, n_sym):
            n_sym = n_sym.title()
            if n_sym in ['C', 'N']:
                if dist < 1.20: return 3.0
                if dist < 1.30: return 2.0
                if dist < 1.36: return 1.5
                if dist < 1.41: return 1.2
                return 1.0
            if n_sym in ['S', 'P']:
                if dist < 1.58: return 2.0
                if dist < 1.66: return 1.5
                return 1.0
            return 1.0

        if sym == 'N':
            total_bo = 0.0
            for n in neighbors:
                n_coord = np.array([n.coordinates.x, n.coordinates.y, n.coordinates.z])
                d = np.linalg.norm(atom_c - n_coord)
                bo = get_n_bond_order(d, n.atomic_symbol)
                total_bo += bo
            
            gap = 3.0 - total_bo
            if gap <= 0.4: needed_h = 0
            elif gap <= 1.4: needed_h = 1
            else: needed_h = 2
            
            needed_h = max(0, min(3, needed_h))
            if n_neigh == 2 and needed_h == 0: return "sp2_carbon", n_neigh
            if n_neigh == 1 and needed_h == 0: return "linear", 1
            target_valency = n_neigh + needed_h
            if n_neigh + needed_h == 4: geom = "tetrahedral"
            elif n_neigh + needed_h == 3: geom = "sp2_carbon"
            else: geom = "linear"
            return geom, target_valency

        elif sym == 'C':
            dists = []
            n_coords = []
            for n in neighbors:
                coord = np.array([n.coordinates.x, n.coordinates.y, n.coordinates.z])
                n_coords.append(coord)
                dists.append(np.linalg.norm(atom_c - coord))
            
            avg_dist = sum(dists)/len(dists) if dists else 0.0
            angles = []
            for i in range(len(n_coords)):
                for j in range(i + 1, len(n_coords)):
                    v1 = n_coords[i] - atom_c
                    v2 = n_coords[j] - atom_c
                    norm = np.linalg.norm(v1) * np.linalg.norm(v2)
                    if norm > 0:
                        dot = np.dot(v1, v2)
                        ang = math.degrees(math.acos(max(-1.0, min(1.0, dot/norm))))
                        angles.append(ang)

            oxy_count = sum(1 for n in neighbors if n.atomic_symbol.strip() == 'O')
            if n_neigh == 1:
                if avg_dist < 1.25: return "linear", 2 
                elif avg_dist < 1.45: return "sp2_carbon", 3 
                else: return "sp3_carbon", 4 
            elif n_neigh == 2:
                if oxy_count >= 1 and avg_dist < 1.40: return "sp2_carbon", 3
                angle = angles[0] if angles else 109.0
                neighbor_syms = [n.atomic_symbol for n in neighbors]
                if 'P' in neighbor_syms:
                     if angle > 165: return "linear", 2
                     else: return "sp3_carbon", 4
                else:
                    if angle > 165: return "linear", 2
                    elif angle > 115: return "sp2_carbon", 3
                    else: return "sp3_carbon", 4
            elif n_neigh == 3:
                if oxy_count >= 2: return "sp2_carbon", 3
                if sum(angles) > 335: return "sp2_carbon", 3
                else: return "sp3_carbon", 4
                
        elif sym == 'S':
            # --- BRANCH 1: Terminal Linear S (e.g., NCS-) ---
            if n_neigh == 1:
                neigh = neighbors[0]
                d = distance(atom, neigh)
                if neigh.atomic_symbol == 'C' and d < 1.75:
                    return "linear", 1 # No hydrogens added
            
            # --- BRANCH 2: Aromatic/Ring Sulfur (e.g., Thiophenes) ---
            # If S has 2 neighbors and the angle is planar (~90-105 deg)
            if n_neigh == 2:
                v1 = get_vec(atom, neighbors[0])
                v2 = get_vec(atom, neighbors[1])
                angle = math.degrees(math.acos(np.dot(normalize(v1), normalize(v2))))
                
                # Thiophene C-S-C angles are typically ~92 degrees
                # If it has 2 heavy neighbors and is part of a ring system, target is 2.
                return "sp2_carbon", 2 

            # --- BRANCH 3: Thiols/Thioethers ---
            return "sp3_carbon", 2 # Target valency of 2
        
        elif sym == 'P':
            # --- BRANCH 1: Trivalent Phosphorus (e.g., Phosphines PR3) ---
            # Target valency is 4 (3 bonds + 1 lone pair slot for H/Metal)
            if n_neigh <= 3:
                return "tetrahedral", 4 
                
            # --- BRANCH 2: Pentavalent Phosphorus (e.g., Phosphates, Phosphonates) ---
            # If already has 4 heavy neighbors (like PO4), valency is satisfied
            if n_neigh == 4:
                return "tetrahedral", 4
                
            return "tetrahedral", 4
        return "tetrahedral", 4
    # --- UNCHANGED: _parse_oxygen_context ---
    def _parse_oxygen_context(self):
        if not self.text_source: return None
        clean_text = self.text_source
        metal_name = ELEMENT_NAMES.get(self.metal_atom.atomic_symbol, "metal").lower()
        parts = clean_text.split(metal_name)
        ligand_sphere_text = parts[0] if len(parts) > 1 else clean_text
        has_aqua = "aqua" in ligand_sphere_text
        has_hydroxo = "hydrox" in ligand_sphere_text
        if has_aqua and not has_hydroxo: return "Water"
        if has_hydroxo and not has_aqua: return "Hydroxide"
        return None
    
    def _calculate_bite_angles(self, connections):
        if len(connections) < 2:
            return []
        metal_coords = np.array([self.metal_atom.coordinates.x, self.metal_atom.coordinates.y, self.metal_atom.coordinates.z])
        vectors = []
        for conn in connections:
            atom = conn.get('atom')
            if not atom or not atom.coordinates:
                continue
            coords = np.array([atom.coordinates.x, atom.coordinates.y, atom.coordinates.z])
            vec = coords - metal_coords
            if np.linalg.norm(vec) > 0:
                vectors.append(vec)
        angles = []
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                v1 = vectors[i]
                v2 = vectors[j]
                norm = np.linalg.norm(v1) * np.linalg.norm(v2)
                if norm > 0:
                    dot = np.dot(v1, v2)
                    ang = math.degrees(math.acos(max(-1.0, min(1.0, dot / norm))))
                    angles.append(float(f"{ang:.3f}"))
        return angles
    
    def _shortest_path_length(self, adjacency, start_idx, end_idx):
        if start_idx == end_idx:
            return 0
        visited = {start_idx}
        queue = [(start_idx, 0)]
        while queue:
            node, dist = queue.pop(0)
            for neighbor in adjacency.get(node, []):
                if neighbor == end_idx:
                    return dist + 1
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, dist + 1))
        return None
    
    def _calculate_chelate_angles(self, connections, component):
        if len(connections) < 2:
            return [], None
        if len(connections) == 2:
            angles = self._calculate_bite_angles(connections)
            terminal_angle = angles[0] if angles else None
            return angles, terminal_angle
        
        adjacency = {}
        for bond in component.bonds:
            a1 = bond.atoms[0].index
            a2 = bond.atoms[1].index
            adjacency.setdefault(a1, []).append(a2)
            adjacency.setdefault(a2, []).append(a1)

        donor_atoms = [c.get('atom') for c in connections if c.get('atom') is not None]
        donor_indices = [a.index for a in donor_atoms if a is not None]
        if len(donor_indices) != len(connections):
            donor_indices = []
        
        terminal_pair = None
        if donor_indices:
            pair_lengths = {}
            for i in range(len(donor_indices)):
                for j in range(i + 1, len(donor_indices)):
                    d1 = donor_indices[i]
                    d2 = donor_indices[j]
                    length = self._shortest_path_length(adjacency, d1, d2)
                    if length is not None:
                        pair_lengths[(d1, d2)] = length
            if pair_lengths:
                terminal_pair = max(pair_lengths, key=pair_lengths.get)
        
        if not terminal_pair:
            donor_indices = sorted([c.get('atom').index for c in connections if c.get('atom') is not None])
            if len(donor_indices) >= 2:
                terminal_pair = (donor_indices[0], donor_indices[-1])

        metal_coords = np.array([self.metal_atom.coordinates.x, self.metal_atom.coordinates.y, self.metal_atom.coordinates.z])
        atom_by_index = {c.get('atom').index: c.get('atom') for c in connections if c.get('atom') is not None}

        def angle_from_indices(i1, i2):
            a1 = atom_by_index.get(i1)
            a2 = atom_by_index.get(i2)
            if not a1 or not a2:
                return None
            v1 = np.array([a1.coordinates.x, a1.coordinates.y, a1.coordinates.z]) - metal_coords
            v2 = np.array([a2.coordinates.x, a2.coordinates.y, a2.coordinates.z]) - metal_coords
            norm = np.linalg.norm(v1) * np.linalg.norm(v2)
            if norm == 0:
                return None
            dot = np.dot(v1, v2)
            ang = math.degrees(math.acos(max(-1.0, min(1.0, dot / norm))))
            return float(f"{ang:.3f}")

        terminal_angle = None
        adjacent_angles = []
        if terminal_pair:
            terminal_angle = angle_from_indices(terminal_pair[0], terminal_pair[1])
            middle_indices = [idx for idx in atom_by_index.keys() if idx not in terminal_pair]
            if len(middle_indices) == 1:
                mid = middle_indices[0]
                for t in terminal_pair:
                    ang = angle_from_indices(t, mid)
                    if ang is not None:
                        adjacent_angles.append(ang)
            else:
                ordered = sorted(atom_by_index.keys())
                for i in range(len(ordered) - 1):
                    ang = angle_from_indices(ordered[i], ordered[i + 1])
                    if ang is not None:
                        adjacent_angles.append(ang)
        
        return adjacent_angles, terminal_angle

    def _process(self):
        mol_copy = self.mol.copy() 
        mol_copy.remove_atoms([a for a in mol_copy.atoms if a.is_metal])
        for component in mol_copy.components:
            if self._is_borate_network(component):
                self._process_borate_component(component)
                continue
            connecting_atoms = []
            for comp_atom in component.atoms:
                orig_atom = self._get_original_atom(comp_atom)
                if not orig_atom: continue
                for neighbor in orig_atom.neighbours:
                    if neighbor.index == self.metal_atom.index:
                        connecting_atoms.append(orig_atom)
                        break
            self._analyze_single_ligand(component, connecting_atoms)

    def _analyze_single_ligand(self, component, connecting_atoms):
        """
        Analyzes the ligand with three possible hydrogen addition paths.
        Path choice: self.addition_method ('None', 'Geometric', 'RDkit')
        """
        non_h_atoms = [a for a in component.atoms if a.atomic_symbol != 'H']
        
        # --- 1. CONNECTIVITY GATHERING ---
        raw_connections = []
        if connecting_atoms:
            center_coords = np.array([self.metal_atom.coordinates.x, self.metal_atom.coordinates.y, self.metal_atom.coordinates.z])
            for catom in connecting_atoms:
                c_coords = np.array([catom.coordinates.x, catom.coordinates.y, catom.coordinates.z])
                raw_connections.append({'atom': catom, 'label': catom.label, 'symbol': catom.atomic_symbol, 'distance': np.linalg.norm(center_coords - c_coords)})

        # Connectivity pruning logic
        indices_to_discard = set()
        n_conn = len(raw_connections)
        for i in range(n_conn):
            for j in range(n_conn):
                if i == j: continue
                if any(b.atoms[0].label in [raw_connections[i]['label'], raw_connections[j]['label']] and 
                       b.atoms[1].label in [raw_connections[i]['label'], raw_connections[j]['label']] for b in component.bonds):
                    if raw_connections[j]['distance'] > (raw_connections[i]['distance'] + self.backbone_tolerance):
                        indices_to_discard.add(j)

        final_connections = [conn for idx, conn in enumerate(raw_connections) if idx not in indices_to_discard]
        dists = [x['distance'] for x in final_connections]
        filtered_connectivity_details = [{'label': x['label'], 'symbol': x['symbol'], 'distance': x['distance']} for x in final_connections]

        is_bridging = any(len([n for n in x['atom'].neighbours if n.is_metal]) > 1 or x['atom'].index in self.known_bridging_indices for x in final_connections)
        if is_bridging: self.has_bridging = True

        denticity = len(final_connections)
        bite_angles = self._calculate_bite_angles(final_connections) if denticity > 1 else []
        chelate_adjacent, chelate_terminal = self._calculate_chelate_angles(final_connections, component) if denticity > 1 else ([], None)

        ligand_dict = {'ligand': '', 'Charge': 0, 'SMILES': '', 'Bridging': is_bridging, 
                       'Connectivity_Details': filtered_connectivity_details, 'Denticity': denticity,
                       'Bite_angle': bite_angles, 'Chelate_angles_adjacent': chelate_adjacent,
                       'Chelate_angle_terminal': chelate_terminal, 'HAD_RDKIT_ISSUE': False}
        local_virtual_hydrogens = []

        # --- BRANCH A: SINGLE ATOM LIGAND (Custom logic for Oxygens remains) ---
        if len(non_h_atoms) == 1 and connecting_atoms:
            atom = non_h_atoms[0]
            orig_atom = self._get_original_atom(atom)
            if orig_atom:
                sym = atom.atomic_symbol
                d = min(dists) if dists else 99.0
                
                if sym == 'O':
                    existing_h = [n for n in orig_atom.neighbours if n.atomic_symbol == 'H']
                    h_count = len(existing_h)
                    species = "Unknown"; needed_h = 0; charge = 0
                    if h_count == 2: species = "Water"; charge = 0
                    elif h_count == 1: species = "Hydroxide"; charge = -1
                    else:
                        metal_sym = self.metal_atom.atomic_symbol
                        oxo_cutoff = actinyl_thresholds.get(metal_sym, 1.98)
                        if d < oxo_cutoff: species = "Oxo"; charge = -2
                        elif is_bridging: species = "Hydroxide"; charge = -1; needed_h = 1
                        else:
                            text_hint = self._parse_oxygen_context()
                            if text_hint == "Water": species = "Water"; charge = 0; needed_h = 2
                            elif text_hint == "Hydroxide": species = "Hydroxide"; charge = -1; needed_h = 1
                            else:
                                if d < 2.25: species = "Hydroxide"; charge = -1; needed_h = 1
                                else: species = "Water"; charge = 0; needed_h = 2

                    ligand_dict.update({'ligand': f"[{species}]", 'Charge': charge, 'SMILES': 'O' if species == 'Oxo' else '[OH-]'})
                    
                    # Logic: Only add H if mode is not 'None'
                    if self.addition_method != 'None' and needed_h > 0:
                        geom = "water_like" if needed_h == 2 else "hydroxyl_like"
                        new_coords = calculate_virtual_hydrogens(orig_atom, [n for n in orig_atom.neighbours], geom, needed_h)
                        for c in new_coords:
                            self.virtual_hydrogens.append(('H', c[0], c[1], c[2]))
                            local_virtual_hydrogens.append(('H', c[0], c[1], c[2], atom.index))

                elif sym in ['F', 'Cl', 'Br', 'I']:
                     ligand_dict.update({'ligand': f"[{sym}-]", 'Charge': -1, 'SMILES': f"[{sym}-]"})

        # --- BRANCH B: POLYATOMIC LIGAND ---
        else:
            # 1. Pre-calculate Hydrogens if using Geometric mode
            if self.addition_method == 'Geometric':
                for atom in component.atoms:
                    orig_atom = self._get_original_atom(atom)
                    if orig_atom and atom.atomic_symbol in ['C', 'N', 'S','P']:
                        # --- NEIGHBOR CHECK ---
                        all_neighbors = [n for n in orig_atom.neighbours if n.coordinates]

                        # --- SPATIAL NEIGHBOR CHECK ---
                        # Instead of just orig_atom.neighbours, we check EVERY atom 
                        # in the cleaned molecule within a 1.3 A radius.
                        #all_neighbors = [
                        #    a for a in self.mol.atoms 
                        #    if a.coordinates and distance(orig_atom, a) < 1.3 and a.index != orig_atom.index
                        #]
                        
                        # Identify heavy neighbors for geometry determination
                        heavy_neighbors = [n for n in all_neighbors if n.atomic_symbol != 'H']
                        
                        # Determine geometry based on the heavy skeleton
                        geometry, target_valency = self._determine_geometry_and_valency(orig_atom, heavy_neighbors)
                        
                        # Calculate needed_h based on physical space occupied
                        needed_h = target_valency - len(all_neighbors)
                    
                        if needed_h > 0:
                            new_coords = calculate_virtual_hydrogens(orig_atom, all_neighbors, geometry, needed_h)
                            
                            for c in new_coords:
                                proposed_xyz = np.array([c[0], c[1], c[2]])
                                is_clash = False
                                
                                # --- FAIL-SAFE: 0.5 A Collision Check ---
                                # Check against all atoms currently in the component
                                for existing_atom in component.atoms:
                                    if existing_atom.coordinates:
                                        ex_xyz = np.array([existing_atom.coordinates.x, 
                                                        existing_atom.coordinates.y, 
                                                        existing_atom.coordinates.z])
                                        dist = np.linalg.norm(proposed_xyz - ex_xyz)
                                        if dist < 0.5:
                                            is_clash = True
                                            #print(f"    [FAIL-SAFE] Blocked H-addition at {dist:.3f} A from {existing_atom.label}")
                                            break
                                
                                if not is_clash:
                                    self.virtual_hydrogens.append(('H', c[0], c[1], c[2]))
                                    local_virtual_hydrogens.append(('H', c[0], c[1], c[2], atom.index))

            # 2. RDKit Processing with Selected Method
            if HAS_RDKIT:
                try:
                    rw_mol = Chem.RWMol()
                    conf = Chem.Conformer()
                    ccdc_to_rdkit_idx = {}

                    # A. Add Heavy Atoms
                    for atom in component.atoms:
                        rd_idx = rw_mol.AddAtom(Chem.Atom(atom.atomic_symbol))
                        ccdc_to_rdkit_idx[atom.index] = rd_idx
                        conf.SetAtomPosition(rd_idx, Geometry.Point3D(atom.coordinates.x, atom.coordinates.y, atom.coordinates.z))
                    # B. Add Hydrogens (Geometric Path)
                    if self.addition_method == 'Geometric':
                        for h_data in local_virtual_hydrogens:
                            h_idx = rw_mol.AddAtom(Chem.Atom('H'))
                            conf.SetAtomPosition(h_idx, Geometry.Point3D(h_data[1], h_data[2], h_data[3]))
                            if len(h_data) > 4:
                                parent_ccdc = h_data[4]
                                if parent_ccdc in ccdc_to_rdkit_idx:
                                    rw_mol.AddBond(ccdc_to_rdkit_idx[parent_ccdc], h_idx, Chem.BondType.SINGLE)

                    rw_mol.AddConformer(conf)
                    # C. Determine Connectivity & Add Hydrogens (RDKit Path)
                    if self.addition_method == 'RDkit':
                        # Must determine connectivity of heavy atoms first
                        rdDetermineBonds.DetermineConnectivity(rw_mol)
                        # AddHs creates coordinate-aware Hydrogens to satisfy valence
                        temp_mol = Chem.AddHs(rw_mol.GetMol(), addCoords=True)
                        rw_mol = Chem.RWMol(temp_mol)
                    else:
                        rdDetermineBonds.DetermineConnectivity(rw_mol)
                    # D. Final Sanitize & SMILES (Hydrogen-aware)
                    base_mol = rw_mol.GetMol()
                    possible_charges = [0, -1, -2, -3, -4, -5, -6, -8, -9, -10, 1, 2]
                    best_mol = None
                    best_energy = float("inf")
                    fallback_mol = None

                    for q in possible_charges:
                        try:
                            temp_mol = Chem.RWMol(base_mol)
                            rdDetermineBonds.DetermineBondOrders(temp_mol, charge=q)
                            candidate_mol = temp_mol.GetMol()
                            Chem.SanitizeMol(candidate_mol)

                            if fallback_mol is None:
                                fallback_mol = Chem.Mol(candidate_mol)

                            energy = float("inf")
                            if candidate_mol.GetNumConformers() > 0:
                                # Fail-safe: MMFF first, then UFF if MMFF cannot run.
                                try:
                                    if AllChem.MMFFHasAllMoleculeParams(candidate_mol):
                                        mmff_props = AllChem.MMFFGetMoleculeProperties(candidate_mol, mmffVariant='MMFF94')
                                        if mmff_props is not None:
                                            mmff_ff = AllChem.MMFFGetMoleculeForceField(candidate_mol, mmff_props)
                                            if mmff_ff is not None:
                                                mmff_ff.Minimize()
                                                energy = mmff_ff.CalcEnergy()
                                except Exception:
                                    pass

                                if not math.isfinite(energy):
                                    try:
                                        if AllChem.UFFHasAllMoleculeParams(candidate_mol):
                                            uff_ff = AllChem.UFFGetMoleculeForceField(candidate_mol)
                                            if uff_ff is not None:
                                                uff_ff.Minimize()
                                                energy = uff_ff.CalcEnergy()
                                    except Exception:
                                        pass
                            if energy < best_energy:
                                best_energy = energy
                                best_mol = Chem.Mol(candidate_mol)

                        except Exception:
                            continue

                    rd_mol = best_mol if best_mol is not None else fallback_mol
                    if rd_mol is None:
                        raise ValueError("No valid charge assignment found for bond-order determination.")

                    Chem.SanitizeMol(rd_mol)
                    ligand_dict.update({'Charge': Chem.GetFormalCharge(rd_mol), 'SMILES': Chem.MolToSmiles(rd_mol, canonical=True)})
                    ligand_dict['ligand'] = ligand_dict['SMILES']

                except Exception:
                    ligand_dict['HAD_RDKIT_ISSUE'] = True
                    ligand_dict['ligand'] = "Polyatomic (Failed)"

        self.ligand_data.append(ligand_dict)

def analyze_ligands(mol, metal_atom, text_source, bridging_indices=None, Hydrogen_Addition_method='Geometric'):
    processor = LigandProcessor(mol, metal_atom, text_source, bridging_indices, 
                                backbone_tolerance=0.3, Hydrogen_Addition_method=Hydrogen_Addition_method)
    
    grouped_data = {}
    for entry in processor.ligand_data:
        name = entry['ligand'] if entry['ligand'] else "Unknown"
        if name not in grouped_data:
            grouped_data[name] = {'ligand': name, 'Charge': entry['Charge'], 'number': 0, 'Bridging': entry['Bridging'], 
                                  'Connectivity_Details': [], 'Denticity_values': [], 'Bite_angle_values': [],
                                  'Chelate_adjacent_values': [], 'Chelate_terminal_values': [],
                                  'HAD_RDKIT_ISSUE': entry.get('HAD_RDKIT_ISSUE', False)}
        group = grouped_data[name]
        group['number'] += 1
        if entry['Bridging']: group['Bridging'] = True
        group['Connectivity_Details'].extend(entry['Connectivity_Details'])
        group['Denticity_values'].append(entry.get('Denticity', 0))
        group['Bite_angle_values'].extend(entry.get('Bite_angle', []))
        group['Chelate_adjacent_values'].extend(entry.get('Chelate_angles_adjacent', []))
        if entry.get('Chelate_angle_terminal') is not None:
            group['Chelate_terminal_values'].append(entry.get('Chelate_angle_terminal'))
    
    final_list = []
    for name, data in grouped_data.items():
        details = data['Connectivity_Details']
        unique_syms = sorted(list(set(d['symbol'] for d in details)))
        avg_per_element = {sym: float(f"{sum([d['distance'] for d in details if d['symbol'] == sym])/len([d['distance'] for d in details if d['symbol'] == sym]):.3f}") for sym in unique_syms}
        avg_denticity = float(f"{sum(data['Denticity_values'])/len(data['Denticity_values']):.3f}") if data['Denticity_values'] else 0.0
        avg_bite_angle = float(f"{sum(data['Bite_angle_values'])/len(data['Bite_angle_values']):.3f}") if data['Bite_angle_values'] else None
        avg_chelate_adjacent = float(f"{sum(data['Chelate_adjacent_values'])/len(data['Chelate_adjacent_values']):.3f}") if data['Chelate_adjacent_values'] else None
        avg_chelate_terminal = float(f"{sum(data['Chelate_terminal_values'])/len(data['Chelate_terminal_values']):.3f}") if data['Chelate_terminal_values'] else None
        final_list.append({
            'ligand': name, 'Charge': data['Charge'], 'number': data['number'], 'Bridging': data['Bridging'],
            'Connecting_atom': [d['label'] for d in details], 'bond_distance': [float(f"{x['distance']:.3f}") for x in details],
            'Connecting_atom_symbol': unique_syms, 'average_bond_distance': avg_per_element,
            'Denticity_values': data['Denticity_values'], 'average_denticity': avg_denticity,
            'Bite_angle_values': data['Bite_angle_values'], 'average_bite_angle': avg_bite_angle,
            'Chelate_angles_adjacent_values': data['Chelate_adjacent_values'], 'average_chelate_adjacent': avg_chelate_adjacent,
            'Chelate_angle_terminal_values': data['Chelate_terminal_values'], 'average_chelate_terminal': avg_chelate_terminal,
            'HAD_RDKIT_ISSUE': data['HAD_RDKIT_ISSUE']
        })
        
    return final_list, processor.has_bridging, processor.oxygen_map, processor.virtual_hydrogens
