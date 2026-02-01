import os
import re

import pandas as pd
import py3Dmol

from .core import process_structure_core
from .ligands import analyze_ligands
from .oxidation import calculate_bvs
from .utils import generate_enriched_xyz_string, view_data2
from .validation import validate_geometry

def parse_manual_input(input_str):
    action = "Keep"
    if "Entry = 'Delete'" in input_str: action = "Delete"
    
    # Extract Deletes (converts to int)
    del_list = []
    del_match = re.search(r"Atom_delete\s*=\s*\[(.*?)\]", input_str)
    if del_match and del_match.group(1).strip():
        del_list = [int(x.strip()) for x in del_match.group(1).split(',')]
        
    # IMPROVED Regex: Finds the content of Atom_average and extracts nested groups
    avg_groups = []
    avg_section = re.search(r"Atom_average\s*=\s*\[(.*)\]", input_str)
    if avg_section:
        content = avg_section.group(1).strip()
        # Find all patterns like [12, 13]
        groups = re.findall(r"\[([\d\s,]+)\]", content)
        for g in groups:
            indices = [int(x.strip()) for x in g.split(',')]
            avg_groups.append(indices)
    return action, {"delete": del_list, "average": avg_groups}

def visualize_manual_disorder(mol, target_metal_label):
    # Initialize a dual-pane viewer
    view = py3Dmol.view(width=1000, height=500, viewergrid=(1,2))
    
    # Generate the base coordinate string
    xyz = generate_enriched_xyz_string(mol, []) 
    
    # --- PANE 1 (LEFT): FULL INDEXING ---
    view.addModel(xyz, 'xyz', viewer=(0,0))
    # Apply standard stick and sphere styling
    view.setStyle({'stick': {'radius':0.1}, 'sphere': {'scale': 0.15}}, viewer=(0,0))
    
    for atom in mol.atoms:
        pos = {'x': atom.coordinates.x, 'y': atom.coordinates.y, 'z': atom.coordinates.z}
        occ = getattr(atom, 'occupancy', 1.0)
        if occ < 0.99:
        # In the left pane, show all labels for identification
            view.addLabel(str(atom.index), {
                'fontSize': 12, 
                'position': pos, 
                'fontColor': 'black', 
                'backgroundColor': 'white', 
                'backgroundOpacity': 0.8,
                'inFront': True,
                'showBackground': False
            }, viewer=(0,0))
        
    # --- PANE 2 (RIGHT): OCCUPANCY-SPECIFIC VIEW ---
    view.addModel(xyz, 'xyz', viewer=(0,1))
    
    # FIX: Ensure bonds/connectivity are visible as in the left image
    view.setStyle({'stick': {'radius':0.1}, 'sphere': {'scale': 0.15}}, viewer=(0,1))
    
    for i, atom in enumerate(mol.atoms):
        occ = getattr(atom, 'occupancy', 1.0)
        
        # Color spheres by occupancy (Blue = Full, Red = Disordered)
        color = 'white' if occ >= 0.99 else 'red'
        view.addStyle({'model': -1, 'serial': i}, {'sphere': {'color': color, 'scale': 0.2}}, viewer=(0,1))
        
        # FIX: Only show labels for atoms with occupancy < 1 in BOTH images
        if occ < 0.99:
            pos = {'x': atom.coordinates.x, 'y': atom.coordinates.y, 'z': atom.coordinates.z}
            label_text = f"{atom.index}: ({occ:.1f})"
            
            # Apply label specifically to disordered atoms in the right pane
            view.addLabel(label_text, {
                'fontSize': 12, 
                'position': pos, 
                'fontColor': 'black',
                'backgroundColor': 'white', 
                'inFront': True,
                'showBackground': False
            }, viewer=(0,1))
            
    view.zoomTo()
    view.show()    

def manual_disorder_html(mol, target_metal_label, width=1000, height=500):
    view = py3Dmol.view(width=width, height=height, viewergrid=(1,2))
    xyz = generate_enriched_xyz_string(mol, [])
    view.addModel(xyz, "xyz", viewer=(0, 0))
    view.setStyle({"stick": {"radius": 0.1}, "sphere": {"scale": 0.15}}, viewer=(0, 0))

    for atom in mol.atoms:
        pos = {"x": atom.coordinates.x, "y": atom.coordinates.y, "z": atom.coordinates.z}
        occ = getattr(atom, "occupancy", 1.0)
        if occ < 0.99:
            view.addLabel(
                str(atom.index),
                {
                    "fontSize": 12,
                    "position": pos,
                    "fontColor": "black",
                    "backgroundColor": "white",
                    "backgroundOpacity": 0.8,
                    "inFront": True,
                    "showBackground": False,
                },
                viewer=(0, 0),
            )

    view.addModel(xyz, "xyz", viewer=(0, 1))
    view.setStyle({"stick": {"radius": 0.1}, "sphere": {"scale": 0.15}}, viewer=(0, 1))

    for i, atom in enumerate(mol.atoms):
        occ = getattr(atom, "occupancy", 1.0)
        color = "white" if occ >= 0.99 else "red"
        view.addStyle(
            {"model": -1, "serial": i},
            {"sphere": {"color": color, "scale": 0.2}},
            viewer=(0, 1),
        )
        if occ < 0.99:
            pos = {"x": atom.coordinates.x, "y": atom.coordinates.y, "z": atom.coordinates.z}
            label_text = f"{atom.index}: ({occ:.1f})"
            view.addLabel(
                label_text,
                {
                    "fontSize": 12,
                    "position": pos,
                    "fontColor": "black",
                    "backgroundColor": "white",
                    "inFront": True,
                    "showBackground": False,
                },
                viewer=(0, 1),
            )

    view.zoomTo()
    return view._make_html()
def manual_refinement_step(TARGET_METAL, num_metal_layers, original_hits, 
                           Hydrogen_Addition_method, Extraction_Cut_off_distances, pkl_name, csv_name,
                           EXTRACTION_METHOD="Topological", Geometric_radius=3.8):
    if not os.path.exists(pkl_name):
        print(f"❌ Could not find result file: {csv_name}")
        return
        
    df = pd.read_pickle(pkl_name)
    to_review = df[(df['VALIDATION_FAILED'] == True) | (df['HAD_RDKIT_ISSUE'] == True)].copy()
    
    if to_review.empty:
        print("✅ No entries found requiring manual correction.")
        return

    for idx, row in to_review.iterrows():
        ref, site = row['CSD ID'], row['Site Label']
        hit = next((h for h in original_hits if h.identifier == ref), None)
        if not hit: continue 

        confirmed = False
        while not confirmed:
            # Phase 1: Visualize current/raw state
            print(f"\n--- [EDITING MODE] {ref} ({site}) ---")
            raw_mol, _, _, _, _, _ = process_structure_core(
                hit.entry, TARGET_METAL, site, num_metal_layers, 
                perform_correction=False,
                EXTRACTION_METHOD=EXTRACTION_METHOD,
                Extraction_Cut_off_distances=Extraction_Cut_off_distances,
                Geometric_radius=Geometric_radius,
            )
            visualize_manual_disorder(raw_mol, site)
            
            user_input = input("Enter Command (e.g., Entry = 'Keep', [Atom_delete= [C1], Atom_average = [[C2, C3]]]) or 'Skip': ")
            
            if user_input.lower() == 'skip': break
            
            action, commands = parse_manual_input(user_input)
            
            if action == 'Delete':
                df = df.drop(idx)
                print(f"🗑️ Entry {ref} deleted from results.")
                confirmed = True
            else:
                # Phase 2: Process with edits and show UPDATE
                print(f"🛠️ Applying edits to {ref}...")
                mol, metal, _, _, bridge, val_f = process_structure_core(
                    hit.entry, TARGET_METAL, site, num_metal_layers, 
                    manual_commands=commands, perform_correction=True,
                    EXTRACTION_METHOD=EXTRACTION_METHOD,
                    Extraction_Cut_off_distances=Extraction_Cut_off_distances,
                    Geometric_radius=Geometric_radius,
                )
                
                # RE-ANALYZE for the preview
                lig_data, has_b, _, v_hs = analyze_ligands(mol, metal, hit.entry.chemical_name, bridge, Hydrogen_Addition_method)
                # RE-VALIDATE including the new hydrogens
                is_actually_valid = validate_geometry(mol, metal, virtual_hydrogens=v_hs, check_rdkit=True)
                val_f = not is_actually_valid # This ensures the preview label is correct

                print(f"Validation Result: {'FAILED' if val_f else 'PASSED'}")
                
                xyz_preview = generate_enriched_xyz_string(mol, v_hs)
                
                print("\n--- [VERIFICATION MODE] Previewing Updated Structure ---")
                print(f"Validation Result: {'FAILED' if val_f else 'PASSED'}")
                view_data2(xyz_preview) # Uses your existing xyz viewer
                
                is_done = input("Are you satisfied with this update? (yes/no): ")
                if is_done.lower() in ['yes', 'y']:
                    # Phase 3: Update DataFrame
                    df.at[idx, 'Coordination Number'] = sum(len(x['Connecting_atom']) for x in lig_data)
                    df.at[idx, 'Ligand_count'] = sum(x['number'] for x in lig_data)
                    df.at[idx, 'AN_n_Ligand_Info'] = str(lig_data)
                    df.at[idx, 'VALIDATION_FAILED'] = val_f
                    df.at[idx, 'HAD_RDKIT_ISSUE'] = any(item.get('HAD_RDKIT_ISSUE', False) for item in lig_data)
                    df.at[idx, 'Bridging_Ligand'] = has_b
                    df.at[idx, 'Total_Ligand_Charge'] = sum(item['Charge'] * item['number'] for item in lig_data)
                    df.at[idx, 'AN_n_XYZ Coordinates'] = xyz_preview
                    print(f"✅ Changes saved for {ref}.")
                    if row['OS_Method'] == "BVS (Geometry)":
                        updated_bvs = calculate_bvs(metal, metal.neighbours)
                        updated_state = f"+{round(updated_bvs)}"
                        df.at[idx, 'Oxidation State'] = updated_state
                    confirmed = True
                else:
                    print("🔄 Returning to editing for the same structure...")

    # Final Save
    df.to_pickle(pkl_name)
    df.to_csv(csv_name, index=False)
    print("\n--- Manual Refinement Completed and Files Updated ---")
