# Function Title: Physical Constants
# Project: Exciton Dynamics
# Author: Doran I. G. Bennett
# Date:January 7, 2016

# Description:
# This is a list of physical constants with the "correct" units.

from numpy import pi


def main():
    pass


if __name__ == "__main__":
    main()

# Calculation Constants
# ---------------------
precision = 1e-8


# Physical Constants
# ------------------
h = 33357.179  # Units: cm^-1*fs
hbar = h / (2 * pi)  # Units: cm^-1*fs/rad
kB = 0.69503  # Units: cm^-1/K
c = 299.79  # Units: nm/fs

# Coupling Constants
dipole_CC = (1 / (4 * pi * 1.5812E-5)) # Units: (cm^-1 *Angstrom^3)/D^2
tresp_CC = 1.1615E5                    # Units: Angstrom*cm^-1)/e^2

# Unit Conversion
# Unit conversion[A][B] converts A<--B
_ev = {"cm-1": 1 / 8065.6, "J": 1 / 1.60218e-19, "meV": 1000.0, "eV": 1.0}
_cm_m1 = {"eV": 8065.6, "J": 8065.6 / 1.60218e-19, "meV": 8.0656, "cm-1": 1}
_j = {
    "cm-1": (1.60218e-19 / 8065.6),
    "eV": (1.60218e-19),
    "meV": ((1.60218e-19) * 1000),
    "J": 1.0,
}
convert_energy_units = {"eV": _ev, "cm-1": _cm_m1, "J": _j}

_eA = {"D": 0.208194}


# Atomic Mass  [units:
# [Copied from BioPython IUPACData.py file]
# Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/ & PyMol
atom_weights = {
    "H": 1.00794,
    "D": 2.01410,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.989770,
    "Mg": 24.3050,
    "Al": 26.981538,
    "Si": 28.0855,
    "P": 30.973761,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955910,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938049,
    "Fe": 55.845,
    "Co": 58.933200,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.39,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.92160,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.80,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.94,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.90545,
    "Ba": 137.327,
    "La": 138.9055,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.24,
    "Pm": 145.0,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92534,
    "Dy": 162.50,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.04,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.9479,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.078,
    "Au": 196.96655,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98038,
    "Po": 208.98,
    "At": 209.99,
    "Rn": 222.02,
    "Fr": 223.02,
    "Ra": 226.03,
    "Ac": 227.03,
    "Th": 232.0381,
    "Pa": 231.03588,
    "U": 238.02891,
    "Np": 237.05,
    "Pu": 244.06,
    "Am": 243.06,
    "Cm": 247.07,
    "Bk": 247.07,
    "Cf": 251.08,
    "Es": 252.08,
    "Fm": 257.10,
    "Md": 258.10,
    "No": 259.10,
    "Lr": 262.11,
    "Rf": 261.11,
    "Db": 262.11,
    "Sg": 266.12,
    "Bh": 264.12,
    "Hs": 269.13,
    "Mt": 268.14,
}

# Atomic Covalent Radii [Units: Angstrom]
# [Copied from BioPython IUPACData.py file]
'''
This is a file transcribed from Table 12 of the journal article below.
The authors note that
* H is from Rowland and Taylor (JPC 1996, 100, 7384) 10.1021/jp953141+
* Be, B, Al, Ca, Ge, Rb, Sr, Sb, Cs, Ba, Bi, Po At, Rn, Fr, Ra are new from the present work
* remainder from Bondi (JPC 1964, 68, 441) 10.1021/j100785a001
'''
dict_covalent_radii = {
 "H": 0.31, "HE": 0.28, "LI": 1.28, "BE": 0.96, "B": 0.85, "C": 0.76,
 "N": 0.71, "O": 0.66, "F": 0.57, "NE": 0.58, "NA": 1.66, "MG": 1.41,
 "AL": 1.21, "SI": 1.11, "P": 1.07, "S": 1.05, "CL": 1.02, "AR": 1.06,
 "K": 2.03, "CA": 1.76, "SC": 1.7, "TI": 1.6, "V": 1.53, "CR": 1.39, "MN": 1.39,
 "FE": 1.32, "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22, "GA": 1.22,
 "GE": 1.2, "AS": 1.19, "SE": 1.2, "BR": 1.2, "KR": 1.16, "RB": 2.2, "SR": 1.95,
 "Y": 1.9, "ZR": 1.75, "NB": 1.64, "MO": 1.54, "TC": 1.47, "RU": 1.46,
 "RH": 1.42, "PD": 1.39, "AG": 1.45, "CD": 1.44, "IN": 1.42, "SN": 1.39,
 "SB": 1.39, "TE": 1.38, "I": 1.39, "XE": 1.4, "CS": 2.44, "BA": 2.15,
 "LA": 2.07, "CE": 2.04, "PR": 2.03, "ND": 2.01, "PM": 1.99, "SM": 1.98,
 "EU": 1.98, "GD": 1.96, "TB": 1.94, "DY": 1.92, "HO": 1.92, "ER": 1.89,
 "TM": 1.9, "YB": 1.87, "LU": 1.87, "HF": 1.75, "TA": 1.7, "W": 1.62,
 "RE": 1.51, "OS": 1.44, "IR": 1.41, "PT": 1.36, "AU": 1.36, "HG": 1.32,
 "TL": 1.45, "PB": 1.46, "BI": 1.48, "PO": 1.4, "AT": 1.5, "RN": 1.5, "FR": 2.6,
 "RA": 2.21, "AC": 2.15, "TH": 2.06, "PA": 2.0, "U": 1.96, "NP": 1.9,
 "PU": 1.87, "AM": 1.8, "CM": 1.69
}

# Atomic Van der Waal Radii [Units: Angstrom]
# Gathered from LANL's Periodic Table
# https://periodic.lanl.gov/list.shtml
dict_vdw_radii = {
    'H': 1.20, 'HE': 1.40,

    'LI': 1.82, 'BE': 1.53, 'B': 1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.35, 'NE': 1.54,

    'NA': 2.27, 'MG': 1.73, 'AL': 1.84, 'SI': 2.10, 'P': 1.80, 'S': 1.80, 'CL': 1.75, 'AR': 1.88,

    'K': 2.75, 'CA': 2.31, 'SC': 2.11, 'TI': 1.87, 'V': 1.79, 'CR': 1.89, 'MN': 1.97, 'FE': 1.94, 'CO': 1.92,
    'NI': 1.63, 'CU': 1.40, 'ZN': 1.39, 'GA': 1.87, 'GE': 2.11, 'AS': 1.85, 'SE': 1.90, 'BR': 1.83, 'KR': 2.02,

    'RB': 3.03, 'SR': 2.49, 'Y': 2.19, 'ZR': 1.86, 'NB': 2.07, 'MO': 2.09, 'TC': 2.09, 'RU': 2.07, 'RH': 1.95,
    'PD': 2.02, 'AG': 1.72, 'CD': 1.58, 'IN': 1.93, 'SN': 2.17, 'SB': 2.06, 'TE': 2.06, 'I': 1.98, 'XE': 2.16,

    'CS': 3.43, 'BA': 2.68, 'LA': 2.40, 'CE': 2.35, 'PR': 2.39, 'ND': 2.29, 'PM': 2.36, 'SM': 2.29, 'EU': 2.33,
    'GD': 2.37, 'TB': 2.21, 'DY': 2.29, 'HO': 2.16, 'ER': 2.35, 'TM': 2.27, 'YB': 2.42, 'LU': 2.21,

    'HF': 2.12, 'TA': 2.17, 'W': 2.10, 'RE': 2.17, 'OS': 2.16, 'IR': 2.02, 'PT': 2.09, 'AU': 1.66, 'HG': 2.09,
    'TL': 1.96, 'PB': 2.02, 'BI': 2.07, 'PO': 1.97, 'AT': 2.02, 'RN': 2.20,

    'FR': 3.48, 'RA': 2.83, 'AC': 2.60, 'TH': 2.37, 'PA': 2.43, 'U': 2.40, 'NP': 2.21, 'PU': 2.43, 'AM': 2.44,
    'CM': 2.45, 'BK': 2.44, 'CF': 2.45, 'ES': 2.45
}

# UFF VDW Parameters
# https://lammpstube.com/wp-content/uploads/2019/10/UFF.pdf
# ------------------
dict_wca_epsilon = {  # Units: cm^-1
    'C': 36.7,
    'O': 21.0,
    'N': 24.1,
    'S': 95.8,
    'MG': 38.8,
    'H': 15.4,
    'FE': 4.6,
    'P': 106.7,
    'NA': 10.5,
    'CA': 83.2,
    'CL': 79.4,
    'MN': 34.5,
}

dict_wca_sigma = {# Units: Angstrom
                    'C': 3.851/2,
                    'O': 3.5/2,
                    'N': 3.66/2,
                    'S': 4.035/2,
                    'MG': 3.021/2, #Assumes Mg2+
                    'H': 2.886/2,
                    'FE': 2.912/2,
                    'P': 4.147/2,
                    'NA': 2.983/2,
                    'CA': 3.399/2, #Assumes Ca2+
                    'CL': 3.947/2
                    }
