import math
import os
import re
import numpy as np

from .constants import ELEMENT_NAMES, BVS_LIBRARY, actinyl_thresholds
from .utils import distance

GOOGLE_API_KEY = ""
try:
    import google.generativeai as genai
    if GOOGLE_API_KEY:
        os.environ['GOOGLE_API_KEY'] = GOOGLE_API_KEY
        genai.configure(api_key=os.environ['GOOGLE_API_KEY'])
        HAS_GEMINI = True
    else:
        HAS_GEMINI = False
except ImportError:
    HAS_GEMINI = False

def resolve_mixed_valence(center_metal, candidate_states, extracted_mol):
    """
    Resolves conflicting oxidation states (e.g. [4, 6], [5, 6], [3, 4]) 
    using Geometry (Actinyl vs Non-Actinyl) and Bond Lengths.
    """
    if not candidate_states: return "?"
    candidate_states = sorted(list(set(candidate_states)))
    if len(candidate_states) == 1: return f"+{candidate_states[0]}"

    # --- 0. DATA PREP ---
    # Approximate midpoints for V vs VI axial bonds for specific actinides.
    # V is generally > threshold, VI is < threshold.
    # derived from general actinide crystal chemistry (U ~1.80, Pu ~1.82-1.84)
  
    # --- 1. GEOMETRY ANALYSIS (Is it Actinyl?) ---
    oxygens = [n for n in center_metal.neighbours if n.atomic_symbol == 'O']
    is_actinyl = False
    max_axial_dist = 0.0
    
    if len(oxygens) >= 2:
        for i in range(len(oxygens)):
            for j in range(i+1, len(oxygens)):
                o1, o2 = oxygens[i], oxygens[j]
                
                v1 = np.array([o1.coordinates.x - center_metal.coordinates.x, o1.coordinates.y - center_metal.coordinates.y, o1.coordinates.z - center_metal.coordinates.z])
                v2 = np.array([o2.coordinates.x - center_metal.coordinates.x, o2.coordinates.y - center_metal.coordinates.y, o2.coordinates.z - center_metal.coordinates.z])
                
                norm1 = np.linalg.norm(v1)
                norm2 = np.linalg.norm(v2)
                
                if norm1 < 0.1 or norm2 < 0.1: continue
                
                dot = np.dot(v1, v2)
                angle = np.degrees(np.arccos(max(-1.0, min(1.0, dot / (norm1 * norm2)))))
                
                # Definition of Actinyl: Linear (~180) and Short (< 2.6 A)
                if angle > 165 and norm1 < 2.6 and norm2 < 2.6:
                    is_actinyl = True
                    # Keep the longer of the two axial bonds for V/VI comparison
                    max_axial_dist = max(norm1, norm2) 
                    break
            if is_actinyl: break

    # --- 2. LOGIC BRANCHING ---

    # Split candidates into "Actinyl-compatible" (V, VI) and "Non-Actinyl" (III, IV)
    actinyl_candidates = [s for s in candidate_states if s in [5, 6]]
    non_actinyl_candidates = [s for s in candidate_states if s in [3, 4]]

    selected_state = None

    if is_actinyl:
        # GEOMETRY SAYS: ACTINYL
        # Logic: If we have actinyl candidates (5 or 6), pick from them.
        # This handles {4, 6} and {4, 5} -> It will ignore the 4 and look at 5 or 6.
        
        if actinyl_candidates:
            if 5 in actinyl_candidates and 6 in actinyl_candidates:
                # SCENARIO: {5, 6}
                # Rule: Longest axial bond is 5.
                # We use the element specific threshold to decide "Long" vs "Short".
                threshold = actinyl_thresholds.get(center_metal.atomic_symbol, 1.82) # Default 1.82 if unknown metal
                
                if max_axial_dist > threshold:
                    selected_state = 5
                else:
                    selected_state = 6
            elif 5 in actinyl_candidates:
                # SCENARIO: {4, 5} -> Geometry confirms Actinyl -> Pick 5
                selected_state = 5
            elif 6 in actinyl_candidates:
                # SCENARIO: {4, 6} -> Geometry confirms Actinyl -> Pick 6
                selected_state = 6
        else:
            # Corner case: Geometry is linear/actinyl, but text only suggested [3, 4].
            # Trust Geometry over Text? usually safe to assume the higher state if geometry forces it.
            # But strictly adhering to text candidates:
            selected_state = max(candidate_states)

    else:
        # GEOMETRY SAYS: NON-ACTINYL (Spherical/Isotropic)
        # Logic: Prefer Non-Actinyl candidates (3, 4).
        # This handles {4, 6} and {4, 5} -> It will ignore the 5/6 and look at 4.
        
        if non_actinyl_candidates:
            if 3 in non_actinyl_candidates and 4 in non_actinyl_candidates:
                # SCENARIO: {3, 4}
                # Rule: Shorter average bond length is 4.
                all_dists = [distance(center_metal, n) for n in center_metal.neighbours]
                avg_dist = sum(all_dists)/len(all_dists) if all_dists else 2.4
                
                # 2.40 is a generic "crossover" point for Pu(III)/Pu(IV)
                if avg_dist < 2.40:
                    selected_state = 4
                else:
                    selected_state = 3
            else:
                # SCENARIO: {4, 6} or {4, 5} but geometry is NOT actinyl -> Pick 4.
                selected_state = non_actinyl_candidates[0] 
        else:
            # Corner case: Geometry is Non-Actinyl, but text only suggested [5, 6].
            # Typically implies Pu(V) without actinyl moiety (rare/unstable) or extraction error.
            # Fallback to lower state.
            selected_state = min(candidate_states)

    return f"+{selected_state}"

def extract_text_oxidation_state(entry, target_metal_symbol):
    text_sources = []
    if entry.chemical_name: text_sources.append(entry.chemical_name)
    if entry.publication and entry.publication.title: text_sources.append(entry.publication.title)
    combined_text = " ".join(text_sources)
    found_states = set()
    roman_map = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7}
    full_name = ELEMENT_NAMES.get(target_metal_symbol, target_metal_symbol)
    targeted_regex = rf"\b({target_metal_symbol}|{full_name})\s*\((I|II|III|IV|V|VI|VII|1\+|2\+|3\+|4\+|5\+|6\+|7\+)\)"
    target_matches = re.findall(targeted_regex, combined_text, re.IGNORECASE)
    targeted_states = set()
    for _, state_str in target_matches: 
        clean = state_str.upper().replace('+', '')
        if clean in roman_map: targeted_states.add(roman_map[clean])
    if targeted_states: return list(targeted_states), combined_text
    prefix_map = {'UNI': 1, 'DI': 2, 'TRI': 3, 'TETRA': 4, 'PENTA': 5, 'HEXA': 6, 'HEPTA': 7, 'OCTA': 8}
    prefix_matches = re.findall(r'\b(uni|di|tri|tetra|penta|hexa|hepta|octa)valent\b', combined_text, re.IGNORECASE)
    valent_states = set()
    for match in prefix_matches:
        if match.upper() in prefix_map: valent_states.add(prefix_map[match.upper()])
    if valent_states: return list(valent_states), combined_text
    matches = re.findall(r'\((I|II|III|IV|V|VI|VII|1\+|2\+|3\+|4\+|5\+|6\+|7\+)\)', combined_text, re.IGNORECASE)
    for match in matches:
        clean = match.upper().replace('+', '')
        if clean in roman_map: found_states.add(roman_map[clean])
    return list(found_states), combined_text

def extract_oxidation_state_with_gemini(abstract, metal_name):
    if not HAS_GEMINI: return None
    model = genai.GenerativeModel('gemini-2.5-flash')
    prompt = f"""Read the following abstract and extract the oxidation state for {metal_name}. Abstract: "{abstract}"
    Rules: 1. Return ONLY the integer value (e.g., +3). 2. If multiple exist, separate by commas. 3. If not mentioned, return "Unknown"."""
    try:
        response = model.generate_content(prompt)
        text = response.text.strip()
        matches = re.findall(r'[+]?(\d)', text)
        if matches: return [int(m) for m in matches]
    except Exception: pass
    return None

def calculate_bvs(metal_atom, neighbors):
    bvs_sum = 0.0
    b_param = 0.37
    metal_sym = metal_atom.atomic_symbol
    metal_params = BVS_LIBRARY.get(metal_sym, {})
    if not metal_params: return 0.0
    for lig in neighbors:
        lig_sym = lig.atomic_symbol
        r0 = metal_params.get(lig_sym)
        c1 = metal_atom.coordinates
        c2 = lig.coordinates
        if c1 is None or c2 is None: continue
        dist = math.sqrt((c1.x - c2.x)**2 + (c1.y - c2.y)**2 + (c1.z - c2.z)**2)
        if r0: bvs_sum += math.exp((r0 - dist) / b_param)
    return bvs_sum

