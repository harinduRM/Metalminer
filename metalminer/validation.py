import math

from .utils import distance

try:
    from rdkit import Chem
    from rdkit import Geometry
    from rdkit.Chem import rdDetermineBonds
    from rdkit import RDLogger
    from rdkit import rdBase
    RDLogger.EnableLog('rdApp.error')
    rdBase.LogToPythonStderr()
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

def validate_with_rdkit(component_atoms, virtual_hydrogens):
    """
    Validates a ligand component by attempting RDKit reconstruction.
    Returns False if sanitization fails (impossible chemistry).
    """
    if not HAS_RDKIT:
        return True # Cannot validate, skip to avoid false negatives
        
    try:
        rw_mol = Chem.RWMol()
        conf = Chem.Conformer()
        ccdc_to_rdkit_idx = {}

        # Add Heavy Atoms
        for atom in component_atoms:
            rd_idx = rw_mol.AddAtom(Chem.Atom(atom.atomic_symbol))
            ccdc_to_rdkit_idx[atom.index] = rd_idx
            conf.SetAtomPosition(rd_idx, Geometry.Point3D(atom.coordinates.x, 
                                                          atom.coordinates.y, 
                                                          atom.coordinates.z))

        # Add Virtual Hydrogens
        for h_data in virtual_hydrogens:
            h_idx = rw_mol.AddAtom(Chem.Atom('H'))
            conf.SetAtomPosition(h_idx, Geometry.Point3D(h_data[1], h_data[2], h_data[3]))
            if len(h_data) > 4:
                parent_ccdc = h_data[4]
                if parent_ccdc in ccdc_to_rdkit_idx:
                    rw_mol.AddBond(ccdc_to_rdkit_idx[parent_ccdc], h_idx, Chem.BondType.SINGLE)

        rw_mol.AddConformer(conf)

        # Determine connectivity and check if it makes chemical sense
        rdDetermineBonds.DetermineConnectivity(rw_mol)
        rd_mol = rw_mol.GetMol()
        # If SanitizeMol fails, it raises an Exception
        Chem.SanitizeMol(rd_mol)
        return True
    
    except Exception:
        return False

def validate_geometry(mol, center_atom, virtual_hydrogens=None, check_rdkit=True, emit=print):
    # 1. Gather all heavy atoms from the molecule
    check_atoms = [a for a in mol.atoms if a.coordinates is not None and distance(center_atom, a) < 20.0]
    
    # 2. Add Virtual Hydrogens to the checklist
    if virtual_hydrogens:
        for vh in virtual_hydrogens:
            # vh is (symbol, x, y, z, parent_index_optional)
            dummy_h = type('DummyAtom', (object,), {
                'atomic_symbol': vh[0],
                'label': 'vH',
                'index': -1,  # Provide a default index for compatibility
                'coordinates': type('Coords', (object,), {'x': vh[1], 'y': vh[2], 'z': vh[3]})
            })
            check_atoms.append(dummy_h)

    n_count = len(check_atoms)   

    # LOOP 1: PAIRWISE CHECKS (Distance)
    for i in range(n_count):
        for j in range(i + 1, n_count):
            a1 = check_atoms[i]
            a2 = check_atoms[j]
            if a1 is a2: continue
            
            # 1. Clean Symbols (Fixes 'O ' vs 'O' error)
            sym1 = a1.atomic_symbol.strip()
            sym2 = a2.atomic_symbol.strip()
            
            d = distance(a1, a2)

            # A. Check for H-H Clashes (Both are Hydrogen)
            if sym1 == 'H' and sym2 == 'H':
                if d < 0.85:  # Tight threshold: H-H shouldn't be this close in coordination sites
                    emit(f"    [VALIDATION FAIL] H-H Overlap: {a1.label}-{a2.label} = {d:.3f} Å")
                    return False

            # B. Check for H-Heavy Clashes (One is Hydrogen)
            elif 'H' in [sym1, sym2]:
                if d < 0.7: # H should not be inside the van der Waals radius of a heavy atom
                    emit(f"    [VALIDATION FAIL] H-Heavy Clash: {a1.label}-{a2.label} = {d:.3f} Å")
                    return False

            # C. Check for Heavy-Heavy Clashes (Neither is Hydrogen)
            else:
                if d < 0.8: 
                    emit(f"    [VALIDATION FAIL] Heavy Atom Clash: {a1.label}-{a2.label} = {d:.2f} Å")
                    return False
                
            # D. O-O Short Contact (Specific check for peroxides)
            if sym1 == 'O' and sym2 == 'O' and d < 1.40:
                emit(f"    [VALIDATION FAIL] Short O-O bond: {d:.3f} Å")
                return False

    # LOOP 2: CARBON ANGLES (Geometry)
    for atom in check_atoms:
        if atom.atomic_symbol.strip() == 'C':
            c_neighbors = []
            has_supporting_neighbor = False
            for potential_n in check_atoms:
                if potential_n.index == atom.index: continue

                d = distance(atom, potential_n)
                n_sym = potential_n.atomic_symbol.strip()
                
                # --- NEW: ISOLATION CHECK ---
                # Check for ANY valid neighbor (C, N, O, H) within 2.0 Å
                if d < 2.0 and n_sym in ['C', 'N', 'O', 'H','P','S','Cl','Br','F']:
                    has_supporting_neighbor = True

                # --- EXISTING: ANGLE GATHERING ---
                # Look for connected carbons for the angle check
                if n_sym == 'C' and 1.1 < d < 1.7: 
                    c_neighbors.append(potential_n)

            # FAIL if the carbon is floating/isolated
            if not has_supporting_neighbor:
                emit(f"    [VALIDATION FAIL] Geometry: Isolated Carbon atom {atom.label} (No C,N,O,H within 2.0Å)")
                return False
            
            if len(c_neighbors) < 2: continue
            
            for i in range(len(c_neighbors)):
                for j in range(i + 1, len(c_neighbors)):
                    n1, n2 = c_neighbors[i], c_neighbors[j]
                    v1 = (n1.coordinates.x - atom.coordinates.x, n1.coordinates.y - atom.coordinates.y, n1.coordinates.z - atom.coordinates.z)
                    v2 = (n2.coordinates.x - atom.coordinates.x, n2.coordinates.y - atom.coordinates.y, n2.coordinates.z - atom.coordinates.z)
                    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
                    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
                    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
                    if mag1 < 0.1 or mag2 < 0.1: continue
                    cos_theta = dot / (mag1 * mag2)
                    cos_theta = max(min(cos_theta, 1.0), -1.0)
                    angle_deg = math.degrees(math.acos(cos_theta))
                    
                    if angle_deg < 90.0: 
                        emit(f"    [VALIDATION FAIL] Geometry: Sharp C-C-C angle ({angle_deg:.1f}°) at {atom.label}")
                        return False

        # --- New RDKit Layer ---
    if check_rdkit and HAS_RDKIT:
        # We must ignore the metal to validate organic ligands
        mol_no_metals = mol.copy()
        mol_no_metals.remove_atoms([a for a in mol_no_metals.atoms if a.is_metal])
        
        for component in mol_no_metals.components:
            # Skip simple inorganic ions/atoms (e.g., Cl-, H2O) 
            # as they rarely fail sanitization unless overlapping
            if len([a for a in component.atoms if a.atomic_symbol != 'H']) < 2:
                continue
            
            # Note: At this stage in process_structure_core, virtual Hs 
            # haven't been calculated yet. We validate the "heavy" skeleton.
            if not validate_with_rdkit(component.atoms, []):
                emit(f"    [VALIDATION FAIL] RDKit could not validate component: {mol.identifier}")
                return False
            
    return True
