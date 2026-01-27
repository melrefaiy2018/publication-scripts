from collections import defaultdict
from copy import deepcopy
import json

def remap_cla_names(data, reverse_mapping, on_collision="keep-both"):
    out = deepcopy(data)

    # 1) Map names (default to original if not found)
    mapped = [reverse_mapping.get(a, a) for a in data["atom"]]

    # 2) Handle collisions
    name_positions = defaultdict(list)
    for i, name in enumerate(mapped):
        name_positions[name].append(i)
    duplicates = {k: v for k, v in name_positions.items() if len(v) > 1}

    if on_collision == "enumerate":
        seen = defaultdict(int)
        mapped_new = []
        for name in mapped:
            seen[name] += 1
            mapped_new.append(name if seen[name] == 1 else f"{name}_{seen[name]}")
        mapped = mapped_new

    elif on_collision == "first-wins":
        keep = set()
        mask = []
        for name in mapped:
            if name not in keep:
                keep.add(name)
                mask.append(True)
            else:
                mask.append(False)
        mapped = [n for n, m in zip(mapped, mask) if m]
        out["q_00"] = [v for v, m in zip(out["q_00"], mask) if m]
        out["q_11"] = [v for v, m in zip(out["q_11"], mask) if m]

    out["atom"] = mapped
    return out, duplicates

def dump_with_wrap(data, file, wrap=8):
    def format_list(lst):
        lines = []
        for i in range(0, len(lst), wrap):
            chunk = lst[i:i+wrap]
            line = ', '.join(json.dumps(x) for x in chunk)
            lines.append(line)
        return '[\n    ' + ',\n    '.join(lines) + '\n  ]'

    def serialize(obj, indent=2):
        if isinstance(obj, dict):
            items = []
            for k, v in obj.items():
                items.append(' ' * indent + json.dumps(k) + ': ' + serialize(v, indent + 2))
            return '{\n' + ',\n'.join(items) + '\n' + ' ' * (indent - 2) + '}'
        elif isinstance(obj, list):
            if len(obj) > wrap and all(isinstance(x, (int, float, str)) for x in obj):
                return format_list(obj)
            else:
                items = []
                for v in obj:
                    items.append(' ' * indent + serialize(v, indent + 2))
                return '[\n' + ',\n'.join(items) + '\n' + ' ' * (indent - 2) + ']'
        else:
            return json.dumps(obj)

    file.write(serialize(data))

# --- Example usage ---
reverse_mapping = {  # (your mapping from non-standard -> standard)
    "C1":"C1","C10":"C10","C11":"C11","C12":"C12","C13":"C13","C14":"C14","C15":"C15","C16":"C16","C17":"C17","C18":"C18","C19":"C19",
    "C1A":"C1A","C1B":"C1B","C1C":"C1C","C1D":"C1D","C2":"C2","C20":"C20","C2A":"C2A","C2B":"C2B","C2C":"C2C","C2D":"C2D","C3":"C3",
    "C3A":"C3A","C3B":"C3B","C3C":"C3C","C3D":"C3D","C4":"C4","C4A":"C4A","C4B":"C4B","C4C":"C4C","C4D":"C4D","C5":"C5","C6":"C6","C7":"C7","C8":"C8","C9":"C9",
    "CAA":"CAA","CAB":"CAB","CAC":"CAC","CAD":"CAD","CBA":"CBA","CBB":"CBB","CBC":"CBC","CBD":"CBD","CED":"CED",
    "CGA":"CGA","CGD":"CGD","CHA":"CHA","CHB":"CHB","CHC":"CHC","CHD":"CHD","CMA":"CMA","CMB":"CMB","CMC":"CMC","CMD":"CMD",
    "H1":"H11","H10":"H61","H11":"H62","H12":"H71","H13":"H72","H14":"H8","H15":"H91","H16":"H92","H17":"H93","H18":"H101","H19":"H102","H2":"H12","H20":"H111",
    "H21":"H112","H22":"H121","H23":"H122","H24":"H13","H25":"H141","H26":"H142","H27":"H143","H28":"H151","H29":"H152","H2A":"H2A","H3":"H2","H30":"H161",
    "H31":"H162","H32":"H171","H33":"H172","H34":"H18","H35":"H191","H36":"H192","H37":"H193","H38":"H202","H39":"H201","H3A":"H3A","H40":"H203",
    "H5":"H41","H6":"H42","H7":"H43","H8":"H51","H9":"H52",
    "HAA1":"HAA1","HAA2":"HAA2","HAB":"HBB","HAC1":"HAC1","HAC2":"HAC2","HBA1":"HBA1","HBA2":"HBA2","HBB1":"HBB1","HBB2":"HBB3",
    "HBC1":"HBC1","HBC2":"HBC2","HBC3":"HBC3","HBD1":"HBD","HED1":"HED1","HED2":"HED2","HED3":"HED3","HHB":"HHB","HHC":"HHC","HHD":"HHD",
    "HMA1":"HMA1","HMA2":"HMA2","HMA3":"HMA3","HMB1":"HMB1","HMB2":"HMB2","HMB3":"HMB3","HMC1":"HMC1","HMC2":"HMC2","HMC3":"HMC3","HMD1":"HMD1","HMD2":"HMD2","HMD3":"HMD3",
    "MG":"MG","NA":"N1A","NB":"N1B","NC":"N1C","ND":"N1D","O1A":"O1A","O1D":"O1D","O2A":"O2A","O2D":"O2D","OBD":"OBD"
}

B3LYP_charges= {
    # Original data
    'MG': (0.912, 0.910), 'CHA': (0.416, 0.434), 'CHB': (-0.386, -0.357),
    'HHB': (0.125, 0.118), 'CHC': (-0.055, -0.127), 'HHC': (0.074, 0.077),
    'CHD': (-0.223, -0.325), 'HHD': (0.147, 0.153), 'NA': (-0.444, -0.452),
    'C1A': (-0.085, -0.060), 'C2A': (-0.318, -0.313), 'H2A': (0.118, 0.118),
    'C3A': (0.427, 0.430), 'H3A': (-0.019, -0.019), 'C4A': (0.300, 0.304),
    'CMA': (-0.584, -0.580), 'HMA1': (0.142, 0.141), 'HMA2': (0.142, 0.141),
    'HMA3': (0.142, 0.141), 'CAA': (0.244, 0.247), 'HAA1': (-0.017, -0.017),
    'HAA2': (-0.017, -0.017), 'CBA': (-0.369, -0.372), 'HBA1': (0.100, 0.101),
    'HBA2': (0.100, 0.101), 'CGA': (0.737, 0.739), 'O1A': (-0.526, -0.526),
    'O2A': (-0.331, -0.331), 'NB': (-0.388, -0.381), 'C1B': (0.124, 0.096),
    'C2B': (0.133, 0.100), 'C3B': (-0.050, -0.078), 'C4B': (0.109, 0.155),
    'CMB': (-0.376, -0.367), 'HMB1': (0.113, 0.110),

    # Extended data
    'HMB2': (0.113, 0.110), 'HMB3': (0.113, 0.110), 'CAB': (0.104, 0.100),
    'HAB': (0.129, 0.128), 'CBB': (-0.358, -0.373), 'HBB1': (0.157, 0.158),
    'HBB2': (0.157, 0.158), 'NC': (-0.347, -0.398), 'C1C': (-0.117, -0.050),
    'C2C': (0.349, 0.356), 'C3C': (-0.253, -0.270), 'C4C': (0.175, 0.288),
    'CMC': (-0.603, -0.603), 'HMC1': (0.164, 0.165), 'HMC2': (0.164, 0.165),
    'HMC3': (0.164, 0.165), 'CAC': (0.145, 0.147), 'HAC1': (0.002, 0.001),
    'HAC2': (0.002, 0.001), 'CBC': (-0.228, -0.228), 'HBC1': (0.061, 0.061),
    'HBC2': (0.061, 0.061), 'HBC3': (0.061, 0.061), 'ND': (-0.383, -0.395),
    'C1D': (0.096, 0.175), 'C2D': (0.237, 0.195), 'C3D': (-0.242, -0.239),
    'C4D': (-0.070, -0.044), 'CMD': (-0.503, -0.494), 'HMD1': (0.147, 0.144),
    'HMD2': (0.147, 0.144), 'HMD3': (0.147, 0.144), 'CAD': (0.673, 0.649),
    'OBD': (-0.477, -0.484), 'CBD': (-0.832, -0.845), 'HBD1': (0.257, 0.261),
    'CGD': (0.827, 0.833), 'O1D': (-0.525, -0.526), 'O2D': (-0.257, -0.258),
    'CED': (-0.272, -0.271), 'HED1': (0.140, 0.140), 'HED2': (0.140, 0.140),
    'HED3': (0.140, 0.140), 'C1': (-0.158, -0.158), 'H1': (0.171, 0.171),
    'H2': (0.171, 0.171)
}

CLA_B3LYP_standardized = {}
CLA_B3LYP_standardized['atom'] = list(B3LYP_charges.keys())
CLA_B3LYP_standardized['q_00'] = [B3LYP_charges[a][0] for a in B3LYP_charges.keys()]
CLA_B3LYP_standardized['q_11'] = [B3LYP_charges[a][1] for a in B3LYP_charges.keys()]

# --- Example usage with your data ---
updated, duplicates = remap_cla_names(CLA_B3LYP_standardized, reverse_mapping, on_collision="keep-both")

# Save to JSON in pretty format with wrapped lists
with open("CLA_B3LYP_standardized.json", "w") as f:
    dump_with_wrap(updated, f, wrap=8)

print("Saved to CLA_B3LYP_standardized.json")

bhhlyp_charges = {
    # Data from the first image
    'MG': (0.964, 0.962), 'CHA': (0.417, 0.430), 'CHB': (-0.432, -0.395),
    'HHB': (0.136, 0.129), 'CHC': (-0.051, -0.122), 'HHC': (0.082, 0.087),
    'CHD': (-0.195, -0.331), 'HHD': (0.153, 0.162), 'NA': (-0.473, -0.484),
    'C1A': (-0.089, -0.078), 'C2A': (-0.305, -0.301), 'H2A': (0.123, 0.123),
    'C3A': (0.417, 0.419), 'H3A': (-0.007, -0.008), 'C4A': (0.331, 0.332),
    'CMA': (-0.615, -0.612), 'HMA1': (0.151, 0.150), 'HMA2': (0.151, 0.150),
    'HMA3': (0.151, 0.150), 'CAA': (0.231, 0.231), 'HAA1': (-0.014, -0.013),
    'HAA2': (-0.014, -0.013), 'CBA': (-0.375, -0.378), 'HBA1': (0.105, 0.106),
    'HBA2': (0.105, 0.106), 'CGA': (0.781, 0.783), 'O1A': (-0.557, -0.558),
    'O2A': (-0.361, -0.362), 'NB': (-0.420, -0.408), 'C1B': (0.162, 0.118),
    'C2B': (0.111, 0.086), 'C3B': (-0.018, -0.054), 'C4B': (0.095, 0.151),
    'CMB': (-0.385, -0.374), 'HMB1': (0.118, 0.115),

    # Newly extracted data
    'HMB2': (0.118, 0.115), 'HMB3': (0.118, 0.115), 'CAB': (-0.136, -0.130),
    'HAB': (0.144, 0.144), 'CBB': (-0.366, -0.382), 'HBB1': (0.168, 0.169),
    'HBB2': (0.168, 0.169), 'NC': (-0.349, -0.417), 'C1C': (-0.121, -0.058),
    'C2C': (0.346, 0.362), 'C3C': (-0.254, -0.271), 'C4C': (0.144, 0.288),
    'CMC': (-0.618, -0.620), 'HMC1': (0.170, 0.171), 'HMC2': (0.170, 0.171),
    'HMC3': (0.170, 0.171), 'CAC': (0.145, 0.146), 'HAC1': (0.004, 0.003),
    'HAC2': (0.004, 0.003), 'CBC': (-0.229, -0.229), 'HBC1': (0.062, 0.062),
    'HBC2': (0.062, 0.062), 'HBC3': (0.062, 0.062), 'ND': (-0.402, -0.426),
    'C1D': (0.076, 0.188), 'C2D': (0.235, 0.186), 'C3D': (-0.240, -0.236),
    'C4D': (-0.072, -0.043), 'CMD': (-0.511, -0.503), 'HMD1': (0.152, 0.150),
    'HMD2': (0.152, 0.150), 'HMD3': (0.152, 0.150), 'CAD': (0.700, 0.686),
    'OBD': (-0.503, -0.506), 'CBD': (-0.844, -0.864), 'HBD1': (0.266, 0.270),
    'CGD': (0.874, 0.881), 'O1D': (-0.556, -0.557), 'O2D': (-0.286, -0.286),
    'CED': (-0.281, -0.281), 'HED1': (0.146, 0.146), 'HED2': (0.146, 0.146),
    'HED3': (0.146, 0.146), 'C1': (-0.162, -0.162), 'H1': (0.180, 0.180),
    'H2': (0.180, 0.180)
}

CLA_BHHLYP_standardized = {}
CLA_BHHLYP_standardized['atom'] = list(bhhlyp_charges.keys())
CLA_BHHLYP_standardized['q_00'] = [bhhlyp_charges[a][0] for a in bhhlyp_charges.keys()]
CLA_BHHLYP_standardized['q_11'] = [bhhlyp_charges[a][1] for a in bhhlyp_charges.keys()]

# --- Example usage with your data ---
updated, duplicates = remap_cla_names(CLA_BHHLYP_standardized, reverse_mapping, on_collision="keep-both")

# Save to JSON in pretty format with wrapped lists
with open("CLA_BHHLYP_standardized.json", "w") as f:
    dump_with_wrap(updated, f, wrap=8)

print("Saved to CLA_BHHLYP_standardized.json")
b65lyp_charges = {
    # Data from the first image
    'MG': (0.994, 0.992), 'CHA': (0.414, 0.429), 'CHB': (-0.464, -0.420),
    'HHB': (0.143, 0.135), 'CHC': (-0.049, -0.116), 'HHC': (0.086, 0.091),
    'CHD': (-0.170, -0.327), 'HHD': (0.155, 0.166), 'NA': (-0.490, -0.498),
    'C1A': (-0.089, -0.084), 'C2A': (-0.301, -0.297), 'H2A': (0.126, 0.127),
    'C3A': (0.414, 0.415), 'H3A': (0.000, -0.001), 'C4A': (0.351, 0.347),
    'CMA': (-0.635, -0.632), 'HMA1': (0.157, 0.156), 'HMA2': (0.157, 0.156),
    'HMA3': (0.157, 0.156), 'CAA': (0.224, 0.222), 
    'HAA1': (-0.011, -0.010),  # Corrected
    'HAA2': (-0.011, -0.010),  # Corrected
    'CBA': (-0.381, -0.384), 'HBA1': (0.108, 0.109), 'HBA2': (0.108, 0.109), 
    'CGA': (0.803, 0.805), 'O1A': (-0.573, -0.573), 'O2A': (-0.375, -0.376), 
    'NB': (-0.439, -0.423), 'C1B': (0.190, 0.137), 'C2B': (0.098, 0.077), 
    'C3B': (0.000, -0.040), 'C4B': (0.086, 0.144), 'CMB': (-0.393, -0.381), 
    'HMB1': (0.121, 0.119),

    # Data from the second image
    'HMB2': (0.121, 0.119), 'HMB3': (0.121, 0.119), 'CAB': (-0.153, -0.147),
    'HAB': (0.152, 0.152), 'CBB': (-0.371, -0.386), 'HBB1': (0.173, 0.175),
    'HBB2': (0.173, 0.175), 'NC': (-0.347, -0.421), 'C1C': (-0.123, -0.065),
    'C2C': (0.345, 0.363), 'C3C': (-0.254, -0.272), 'C4C': (0.122, 0.280),
    'CMC': (-0.629, -0.631), 'HMC1': (0.173, 0.174), 'HMC2': (0.173, 0.174),
    'HMC3': (0.173, 0.174), 'CAC': (0.143, 0.145), 'HAC1': (0.006, 0.005),
    'HAC2': (0.006, 0.005), 'CBC': (-0.233, -0.233), 'HBC1': (0.063, 0.063),
    'HBC2': (0.063, 0.063), 'HBC3': (0.063, 0.063), 'ND': (-0.409, -0.437),
    'C1D': (0.058, 0.187), 'C2D': (0.238, 0.187), 'C3D': (-0.241, -0.238),
    'C4D': (-0.071, -0.042), 'CMD': (-0.520, -0.512), 'HMD1': (0.156, 0.153),
    'HMD2': (0.156, 0.153), 'HMD3': (0.156, 0.153), 'CAD': (0.712, 0.702),
    'OBD': (-0.517, -0.519), 'CBD': (-0.852, -0.873), 'HBD1': (0.270, 0.275),
    'CGD': (0.897, 0.904), 'O1D': (-0.571, -0.573), 'O2D': (-0.299, -0.299),
    'CED': (-0.288, -0.289), 'HED1': (0.151, 0.151), 'HED2': (0.151, 0.151),
    'HED3': (0.151, 0.151), 'C1': (-0.168, -0.168), 'H1': (0.184, 0.184),
    'H2': (0.185, 0.185)
}
CLA_B65LYP_standardized = {}
CLA_B65LYP_standardized['atom'] = list(b65lyp_charges.keys())
CLA_B65LYP_standardized['q_00'] = [b65lyp_charges[a][0] for a in b65lyp_charges.keys()]
CLA_B65LYP_standardized['q_11'] = [b65lyp_charges[a][1] for a in b65lyp_charges.keys()]



# --- Example usage with your data ---
updated, duplicates = remap_cla_names(CLA_B65LYP_standardized, reverse_mapping, on_collision="keep-both")

# Save to JSON in pretty format with wrapped lists
with open("CLA_B65LYP_standardized.json", "w") as f:
    dump_with_wrap(updated, f, wrap=8)

print("Saved to CLA_B65LYP_standardized.json")
