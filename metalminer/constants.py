"""Constants for the MetalMiner project."""

ELEMENT_NAMES = {
    # --- Lanthanides (La-Lu) ---
    'La': 'Lanthanum', 'Ce': 'Cerium', 'Pr': 'Praseodymium', 'Nd': 'Neodymium', 'Pm': 'Promethium', 
    'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium', 'Tb': 'Terbium', 'Dy': 'Dysprosium', 
    'Ho': 'Holmium', 'Er': 'Erbium', 'Tm': 'Thulium', 'Yb': 'Ytterbium', 'Lu': 'Lutetium',
    
    # --- Actinides (Ac-Lr) ---
    'Ac': 'Actinium', 'Th': 'Thorium', 'Pa': 'Protactinium', 'U': 'Uranium', 'Np': 'Neptunium', 
    'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium', 'Cf': 'Californium', 
    'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium',
    
    # --- Transition Metals: Period 4 ---
    'Sc': 'Scandium',  'Ti': 'Titanium',  'V': 'Vanadium',   'Cr': 'Chromium',  'Mn': 'Manganese',
    'Fe': 'Iron',      'Co': 'Cobalt',    'Ni': 'Nickel',    'Cu': 'Copper',    'Zn': 'Zinc',
    
    # --- Transition Metals: Period 5 ---
    'Y': 'Yttrium',    'Zr': 'Zirconium', 'Nb': 'Niobium',   'Mo': 'Molybdenum','Tc': 'Technetium',
    'Ru': 'Ruthenium', 'Rh': 'Rhodium',   'Pd': 'Palladium', 'Ag': 'Silver',    'Cd': 'Cadmium',
    
    # --- Transition Metals: Period 6 ---
    'Hf': 'Hafnium',   'Ta': 'Tantalum',  'W': 'Tungsten',   'Re': 'Rhenium',   'Os': 'Osmium',
    'Ir': 'Iridium',   'Pt': 'Platinum',  'Au': 'Gold',      'Hg': 'Mercury',
    
    # --- Transition Metals: Period 7 (Superheavy) ---
    'Rf': 'Rutherfordium', 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium',
    'Mt': 'Meitnerium',    'Ds': 'Darmstadtium', 'Rg': 'Roentgenium', 'Cn': 'Copernicium'
}

# BVS data is adopted from https://doi.org/10.1107/S0108768190011041
BVS_LIBRARY = {
    # --- Lanthanides ---
    'La': {'O': 2.172, 'N': 2.22, 'F': 2.08, 'Cl': 2.60, 'S': 2.57, 'C': 2.30},
    'Ce': {'O': 2.151, 'N': 2.20, 'F': 2.06, 'Cl': 2.58, 'S': 2.55, 'C': 2.28},
    'Pr': {'O': 2.138, 'N': 2.18, 'F': 2.05, 'Cl': 2.56, 'S': 2.54, 'C': 2.26},
    'Nd': {'O': 2.117, 'N': 2.16, 'F': 2.03, 'Cl': 2.55, 'S': 2.53, 'C': 2.25},
    'Pm': {'O': 2.100, 'N': 2.15, 'F': 2.02, 'Cl': 2.54, 'S': 2.51, 'C': 2.23},
    'Sm': {'O': 2.088, 'N': 2.14, 'F': 2.00, 'Cl': 2.53, 'S': 2.50, 'C': 2.21},
    'Eu': {'O': 2.064, 'N': 2.13, 'F': 1.98, 'Cl': 2.52, 'S': 2.49, 'C': 2.15},
    'Gd': {'O': 2.065, 'N': 2.11, 'F': 1.96, 'Cl': 2.50, 'S': 2.47, 'C': 2.15},
    'Tb': {'O': 2.048, 'N': 2.09, 'F': 1.94, 'Cl': 2.49, 'S': 2.45, 'C': 2.13},
    'Dy': {'O': 2.033, 'N': 2.07, 'F': 1.92, 'Cl': 2.47, 'S': 2.43, 'C': 2.11},
    'Ho': {'O': 2.022, 'N': 2.06, 'F': 1.91, 'Cl': 2.46, 'S': 2.42, 'C': 2.10},
    'Er': {'O': 2.005, 'N': 2.04, 'F': 1.89, 'Cl': 2.45, 'S': 2.41, 'C': 2.08},
    'Tm': {'O': 1.991, 'N': 2.03, 'F': 1.88, 'Cl': 2.44, 'S': 2.40, 'C': 2.07},
    'Yb': {'O': 1.974, 'N': 2.02, 'F': 1.87, 'Cl': 2.43, 'S': 2.38, 'C': 2.05},
    'Lu': {'O': 1.961, 'N': 2.00, 'F': 1.85, 'Cl': 2.41, 'S': 2.37, 'C': 2.04},

    # --- Actinides ---
    'Ac': {'O': 2.24, 'N': 2.30, 'F': 2.15, 'Cl': 2.65, 'S': 2.65, 'C': 2.35},
    'Th': {'O': 2.18, 'N': 2.25, 'F': 2.10, 'Cl': 2.60, 'S': 2.60, 'C': 2.30},
    'Pa': {'O': 2.12, 'N': 2.20, 'F': 2.06, 'Cl': 2.55, 'S': 2.55, 'C': 2.25},
    'U':  {'O': 2.059, 'N': 2.12, 'F': 1.98, 'Cl': 2.43, 'S': 2.43, 'C': 2.14},
    'Np': {'O': 2.05, 'N': 2.12, 'F': 2.00, 'Cl': 2.48, 'S': 2.48, 'C': 2.18},
    'Pu': {'O': 2.03, 'N': 2.10, 'F': 1.98, 'Cl': 2.46, 'S': 2.46, 'C': 2.16},
    'Am': {'O': 2.02, 'N': 2.08, 'F': 1.96, 'Cl': 2.45, 'S': 2.45, 'C': 2.15},
    'Cm': {'O': 2.00, 'N': 2.06, 'F': 1.95, 'Cl': 2.44, 'S': 2.44, 'C': 2.14},
    'Bk': {'O': 2.06, 'N': 2.10, 'F': 1.98, 'Cl': 2.45, 'S': 2.45, 'C': 2.15},
    'Cf': {'O': 1.96, 'N': 2.02, 'F': 1.92, 'Cl': 2.40, 'S': 2.40, 'C': 2.12},
    'Es': {'O': 1.95, 'N': 2.00, 'F': 1.90, 'Cl': 2.39, 'S': 2.39, 'C': 2.10},
    'Fm': {'O': 1.94, 'N': 1.98, 'F': 1.88, 'Cl': 2.38, 'S': 2.38, 'C': 2.08},
    
    # --- Transition Metals (Period 4) ---
    # Values mostly for M(II) or M(III) where typical
    'Sc': {'O': 1.849, 'N': 1.88, 'F': 1.79, 'Cl': 2.26, 'S': 2.28, 'C': 2.05}, # Sc(III)
    'Ti': {'O': 1.815, 'N': 1.86, 'F': 1.77, 'Cl': 2.22, 'S': 2.23, 'C': 1.96}, # Ti(IV)
    'V':  {'O': 1.803, 'N': 1.84, 'F': 1.73, 'Cl': 2.22, 'S': 2.21, 'C': 1.95}, # V(V)
    'Cr': {'O': 1.760, 'N': 1.78, 'F': 1.71, 'Cl': 2.13, 'S': 2.16, 'C': 1.88}, # Cr(III)
    'Mn': {'O': 1.790, 'N': 1.81, 'F': 1.72, 'Cl': 2.15, 'S': 2.20, 'C': 1.95}, # Mn(II)
    'Fe': {'O': 1.759, 'N': 1.77, 'F': 1.68, 'Cl': 2.10, 'S': 2.15, 'C': 1.70}, # Fe(III)
    'Co': {'O': 1.692, 'N': 1.73, 'F': 1.65, 'Cl': 2.05, 'S': 2.09, 'C': 1.75}, # Co(II)
    'Ni': {'O': 1.654, 'N': 1.69, 'F': 1.62, 'Cl': 2.02, 'S': 2.06, 'C': 1.78}, # Ni(II)
    'Cu': {'O': 1.679, 'N': 1.71, 'F': 1.58, 'Cl': 2.00, 'S': 2.05, 'C': 1.65}, # Cu(II)
    'Zn': {'O': 1.704, 'N': 1.75, 'F': 1.64, 'Cl': 2.04, 'S': 2.09, 'C': 1.80}, # Zn(II)

    # --- Transition Metals (Period 5) ---
    'Y':  {'O': 2.019, 'N': 2.07, 'F': 1.93, 'Cl': 2.45, 'S': 2.45, 'C': 2.18}, # Y(III)
    'Zr': {'O': 1.928, 'N': 1.98, 'F': 1.87, 'Cl': 2.33, 'S': 2.34, 'C': 2.07}, # Zr(IV)
    'Nb': {'O': 1.911, 'N': 1.96, 'F': 1.84, 'Cl': 2.29, 'S': 2.31, 'C': 2.05}, # Nb(V)
    'Mo': {'O': 1.907, 'N': 1.92, 'F': 1.83, 'Cl': 2.26, 'S': 2.28, 'C': 2.00}, # Mo(VI)
    'Tc': {'O': 1.874, 'N': 1.90, 'F': 1.80, 'Cl': 2.22, 'S': 2.22, 'C': 2.00}, # Tc(VII) est
    'Ru': {'O': 1.770, 'N': 1.82, 'F': 1.76, 'Cl': 2.16, 'S': 2.17, 'C': 1.90}, # Ru(IV)
    'Rh': {'O': 1.730, 'N': 1.78, 'F': 1.73, 'Cl': 2.11, 'S': 2.12, 'C': 1.85}, # Rh(III)
    'Pd': {'O': 1.650, 'N': 1.69, 'F': 1.67, 'Cl': 2.06, 'S': 2.07, 'C': 1.80}, # Pd(II)
    'Ag': {'O': 1.842, 'N': 1.88, 'F': 1.77, 'Cl': 2.17, 'S': 2.19, 'C': 1.98}, # Ag(I)
    'Cd': {'O': 1.904, 'N': 1.94, 'F': 1.82, 'Cl': 2.29, 'S': 2.32, 'C': 2.05}, # Cd(II)

    # --- Transition Metals (Period 6) ---
    'Hf': {'O': 1.910, 'N': 1.97, 'F': 1.86, 'Cl': 2.32, 'S': 2.33, 'C': 2.06}, # Hf(IV)
    'Ta': {'O': 1.920, 'N': 1.97, 'F': 1.85, 'Cl': 2.31, 'S': 2.32, 'C': 2.05}, # Ta(V)
    'W':  {'O': 1.917, 'N': 1.95, 'F': 1.83, 'Cl': 2.27, 'S': 2.29, 'C': 2.02}, # W(VI)
    'Re': {'O': 1.900, 'N': 1.94, 'F': 1.83, 'Cl': 2.23, 'S': 2.25, 'C': 2.00}, # Re(VII)
    'Os': {'O': 1.780, 'N': 1.83, 'F': 1.79, 'Cl': 2.19, 'S': 2.20, 'C': 1.95}, # Os(IV)
    'Ir': {'O': 1.750, 'N': 1.80, 'F': 1.77, 'Cl': 2.16, 'S': 2.17, 'C': 1.90}, # Ir(IV)
    'Pt': {'O': 1.660, 'N': 1.71, 'F': 1.67, 'Cl': 2.08, 'S': 2.09, 'C': 1.85}, # Pt(II)
    'Au': {'O': 1.690, 'N': 1.70, 'F': 1.65, 'Cl': 2.05, 'S': 2.06, 'C': 1.70}, # Au(III)
    'Hg': {'O': 1.970, 'N': 2.00, 'F': 1.90, 'Cl': 2.34, 'S': 2.37, 'C': 2.05}, # Hg(II)
}


actinyl_thresholds = {'U': 1.98, 'Np': 1.98, 'Pu': 1.98, 'Am': 1.98,'Cm': 1.98}
#=============================================================================

__all__ = ['ELEMENT_NAMES', 'BVS_LIBRARY', 'actinyl_thresholds']

