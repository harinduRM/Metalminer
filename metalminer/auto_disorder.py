import numpy as np
from ccdc import molecule
from scipy.spatial import cKDTree

def resolve_cluster_fast(cluster_data, method="Hybrid"):
    """
    Resolves a specific cluster of atoms based on the selected method.
    Methods: 'Hybrid', 'Major_componet', 'Minor_componet', 'Average'
    """
    if not cluster_data: return None, [], None
    if len(cluster_data) == 1: return cluster_data[0][0], [], None

    # Sort descending: Index 0 is Best (Major), Index -1 is Worst (Minor)
    sorted_data = sorted(cluster_data, key=lambda x: x[1], reverse=True)
    
    # Identify Major (A) and Runner-up (B)
    major_atom = sorted_data[0][0]
    major_occ = sorted_data[0][1]
    major_coords = sorted_data[0][2]
    
    second_occ = sorted_data[1][1]
    second_coords = sorted_data[1][2]

    # Check for 50/50 split (approximate equality)
    occ_diff = abs(major_occ - second_occ)
    is_split_site = occ_diff < 0.05

    winner_atom = major_atom
    atoms_to_remove = []
    new_coords = None

    # --- LOGIC BRANCHING ---

    if method == "Hybrid":
        # Default behavior: Merge if 50/50, otherwise keep Major
        atoms_to_remove = [x[0] for x in sorted_data[1:]]
        if is_split_site:
            avg = (major_coords + second_coords) / 2.0
            new_coords = molecule.Coordinates(avg[0], avg[1], avg[2])

    elif method == "Major_componet":
        # Always keep Major. Never move.
        atoms_to_remove = [x[0] for x in sorted_data[1:]]
        # winner_atom is already major_atom, new_coords remains None

    elif method == "Minor_componet":
        if is_split_site:
            # 0.50 vs 0.50 -> Keep A (Major)
            winner_atom = major_atom
            atoms_to_remove = [x[0] for x in sorted_data[1:]]
        else:
            # 0.80 vs 0.20 -> Keep B (Minor)
            # The minor component is the last item in the sorted list
            minor_data = sorted_data[-1]
            winner_atom = minor_data[0]
            # Remove everything else (including the Major one)
            atoms_to_remove = [x[0] for x in sorted_data[:-1]]

    elif method == "Average":
        # Weighted Average
        atoms_to_remove = [x[0] for x in sorted_data[1:]]
        
        # Calculate weighted center of mass for the whole cluster
        total_occ = sum(x[1] for x in sorted_data)
        if total_occ > 0:
            # FIX: Use index [0], [1], [2] instead of .x, .y, .z
            w_x = sum(x[2][0] * x[1] for x in sorted_data) / total_occ
            w_y = sum(x[2][1] * x[1] for x in sorted_data) / total_occ
            w_z = sum(x[2][2] * x[1] for x in sorted_data) / total_occ
            new_coords = molecule.Coordinates(w_x, w_y, w_z)

    return winner_atom, atoms_to_remove, new_coords

def clean_small_cluster(mol, atom_to_layer_map, perform_correction=True, disorder_method="Hybrid"):
    atoms = mol.atoms
    n_atoms = len(atoms)
    
    # 1. Bulk Extract
    coords = np.array([a.coordinates for a in atoms])
    occupancies = np.array([getattr(a, 'occupancy', 1.0) for a in atoms])
    symbols = [a.atomic_symbol for a in atoms]
    is_metals = [a.is_metal for a in atoms] 
    
    # 2. Global Overlap Search
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=2.2) 
    
    if not pairs: return mol, False
    
    # 3. Build Clusters
    adj = [[] for _ in range(n_atoms)]
    active_indices = set()
    found_disorder = False
    
    for i, j in pairs:
        if symbols[i] != symbols[j]: continue
        dist = np.linalg.norm(coords[i] - coords[j])
        is_disorder_pair = False
        
        # RULE A: Impossible Overlap
        if dist < 0.9:
            is_disorder_pair = True
        # RULE B: Heavy Metal Ghost
        elif is_metals[i] and dist < 2.0:
            if occupancies[i] < 1.0 and occupancies[j] < 1.0:
                 is_disorder_pair = True

        if is_disorder_pair:
            adj[i].append(j)
            adj[j].append(i)
            active_indices.add(i)
            active_indices.add(j)
            found_disorder = True

    # --- RULE C: PARTIAL OCCUPANCY CHECK (New) ---
    # Any atom with occupancy < 1.0 is flagged as disorder.
    # We add it to active_indices so it gets passed to the resolver 
    # (even if it has no neighbors, the resolver will just keep it, but the flag is set).
    for i in range(n_atoms):
        if occupancies[i] < 1.0:
            found_disorder = True
            active_indices.add(i)

    if not found_disorder: return mol, False 
    
    if not perform_correction:
        return mol, True 
    
# 4. Resolve
    visited = [False] * n_atoms
    atoms_to_delete = []
    coord_updates = {} 
    
    for i in active_indices:
        if visited[i]: continue
        component_indices = [i]
        queue = [i]
        visited[i] = True
        while queue:
            curr = queue.pop(0)
            for neighbor in adj[curr]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    component_indices.append(neighbor)
                    queue.append(neighbor)
        
        cluster_data = []
        for idx in component_indices:
            cluster_data.append((atoms[idx], occupancies[idx], coords[idx]))
            
        # --- PASS THE METHOD HERE ---
        winner, losers, new_coords = resolve_cluster_fast(cluster_data, method=disorder_method)
        
        if new_coords:
            coord_updates[winner] = new_coords
        atoms_to_delete.extend(losers)

    # 5. Apply Updates
    for atom_obj, new_xyz in coord_updates.items():
        atom_obj.coordinates = new_xyz
        
    if atoms_to_delete:
        mol.remove_atoms(atoms_to_delete)
        
    return mol, True

