import math
import numpy as np
import py3Dmol
import requests

def distance(a1, a2):
    """Calculates Euclidean distance between two CCDC atoms."""
    return math.sqrt(
        (a1.coordinates.x - a2.coordinates.x)**2 + 
        (a1.coordinates.y - a2.coordinates.y)**2 + 
        (a1.coordinates.z - a2.coordinates.z)**2
    )

def get_vec(a1, a2):
    """Return vector from a1 to a2."""
    if hasattr(a1, 'coordinates'): c1 = np.array([a1.coordinates.x, a1.coordinates.y, a1.coordinates.z])
    else: c1 = np.array(a1)
    if hasattr(a2, 'coordinates'): c2 = np.array([a2.coordinates.x, a2.coordinates.y, a2.coordinates.z])
    else: c2 = np.array(a2)
    return c2 - c1

def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-6 else np.array([0, 0, 1])

def view_data2(AN_n_xyz_block):
    view = py3Dmol.view(width=400, height=400)
    view.addModel(AN_n_xyz_block, 'xyz',)
    view.setStyle({'stick': {}, 'sphere': {'scale': 0.2}}) 

    view.addPropertyLabels("atom", {}, {
        'fontColor': 'black',
        'fontOpacity': 1.0, 
        'backgroundColor': 'white', 
        'backgroundOpacity': 0.5, 
        'fontSize': 10
    })

    view.zoomTo()
    view.show()

def make_xyz_viewer_html(AN_n_xyz_block, width=400, height=400):
    view = py3Dmol.view(width=width, height=height)
    view.addModel(AN_n_xyz_block, "xyz")
    view.setStyle({"stick": {}, "sphere": {"scale": 0.2}})
    view.addPropertyLabels(
        "atom",
        {},
        {
            "fontColor": "black",
            "fontOpacity": 1.0,
            "backgroundColor": "white",
            "backgroundOpacity": 0.5,
            "fontSize": 10,
        },
    )
    view.zoomTo()
    return view._make_html()

def generate_enriched_xyz_string(mol, virtual_hydrogens):
    """
    Generates an XYZ block including original atoms plus virtual hydrogens.
    
    FIX: Deduplicates atoms based on coordinates and symbol to prevent 
    'ghost' or duplicate atoms if the extraction process added them twice.
    """
    valid_atoms = []
    seen_sites = set()
    
    # 1. Deduplicate Real Atoms
    for atom in mol.atoms:
        if atom.coordinates is None: continue
        
        # Create a signature based on Symbol and Position (rounded to 3 decimals)
        # This handles cases where 'I2' is added twice at the same spot.
        site_sig = (
            atom.atomic_symbol, 
            round(atom.coordinates.x, 3), 
            round(atom.coordinates.y, 3), 
            round(atom.coordinates.z, 3)
        )
        
        if site_sig not in seen_sites:
            seen_sites.add(site_sig)
            valid_atoms.append(atom)
            
    total_count = len(valid_atoms) + len(virtual_hydrogens)
    
    lines = [str(total_count), f"{mol.identifier} - Enriched"]
    
    # 2. Print Real Atoms
    for atom in valid_atoms:
        lines.append(f"{atom.atomic_symbol:3s} {atom.coordinates.x:10.5f} {atom.coordinates.y:10.5f} {atom.coordinates.z:10.5f}")
    
    # 3. Print Virtual Hydrogens
    for h in virtual_hydrogens:
        lines.append(f"{h[0]:3s} {h[1]:10.5f} {h[2]:10.5f} {h[3]:10.5f}")
        
    return "\n".join(lines)

def get_abstract_from_openalex(doi):
    url = f"https://api.openalex.org/works/https://doi.org/{doi}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code != 200: return None
        data = response.json()
        index = data.get('abstract_inverted_index')
        if index:
            abstract_map = {}
            for word, positions in index.items():
                for pos in positions: abstract_map[pos] = word
            return " ".join([abstract_map[i] for i in sorted(abstract_map.keys())])
    except Exception: pass
    return None

def get_nuclearity_label(count):
    if count == 1: return "Monomer"
    if count == 2: return "Dimer"
    if count == 3: return "Trimer"
    if count == 4: return "Tetramer"
    return "Cluster/Oligomer"
