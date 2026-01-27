# ======================
# Mohamed, June 27,2024:
# ======================

import numpy as np
from Bio.PDB import PDBParser, Superimposer

def align_molecules(pdb_file1, pdb_file2):
    """
    Aligns two molecules (PDB files) based on their atoms' positions.
    
    Parameters:
    - pdb_file1 (str): Path to the first PDB file.
    - pdb_file2 (str): Path to the second PDB file.
    
    Returns:
    - aligned_atoms1 (list): List of aligned Atom objects from molecule 1.
    - aligned_atoms2 (list): List of aligned Atom objects from molecule 2.
    """
    # Load atoms from PDB files
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure("molecule1", pdb_file1)
    structure2 = parser.get_structure("molecule2", pdb_file2)
    
    # Extract atoms from structures
    atoms1 = []
    atoms2 = []
    for model in structure1:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms1.append(atom)
    
    for model in structure2:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms2.append(atom)
    
    # Align structures using Superimposer
    super_imposer = Superimposer()
    # Extract coordinate sets for superimposition
    coordset1 = np.array([atom.get_coord() for atom in atoms1])
    coordset2 = np.array([atom.get_coord() for atom in atoms2])
    
    # Perform superimposition
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(atoms2)
    
    # Return aligned atoms
    aligned_atoms1 = atoms1
    aligned_atoms2 = atoms2
    
    return aligned_atoms1, aligned_atoms2


def load_pdb(file_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("molecule", file_path)
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom)
    return atoms

def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1.get_coord() - atom2.get_coord())

def identify_heavy_atoms(atoms):
    heavy_atoms = []
    for atom in atoms:
        if atom.element != 'H':
            heavy_atoms.append(atom)
    return heavy_atoms

def identify_hydrogen_atoms(heavy_atom):
    hydrogen_atoms = []
    for neighbor in heavy_atom.get_neighbors():
        if neighbor.element == 'H':
            hydrogen_atoms.append(neighbor)
    return hydrogen_atoms

def map_atoms(atoms1, atoms2):
    mapping = {}
    for atom1 in atoms1:
        min_distance = float('inf')
        closest_atom = None
        for atom2 in atoms2:
            distance = calculate_distance(atom1, atom2)
            if distance < min_distance:
                min_distance = distance
                closest_atom = atom2
        mapping[atom1.get_id()] = closest_atom.get_id()
    return mapping

def are_atoms_bonded(atom1, atom2, threshold_distance=2.0):
    # Example: Check if the distance between atom1 and atom2 is below a threshold
    distance = calculate_distance(atom1, atom2)  # Implement calculate_distance as per your requirements
    return distance <= threshold_distance


def identify_hydrogen_atoms(heavy_atom, atoms):
    hydrogen_atoms = []
    for atom in atoms:
        # Example: Check if the atom is a hydrogen atom and bonded to the heavy_atom
        if atom.element == 'H' and are_atoms_bonded(heavy_atom, atom):
            hydrogen_atoms.append(atom)
    return hydrogen_atoms



def calculate_distance(atom1, atom2):
    pos1 = np.array(atom1.get_coord())  
    pos2 = np.array(atom2.get_coord())
    return np.linalg.norm(pos1 - pos2)

def map_hydrogen_atoms(heavy_atoms1, heavy_atoms2, atoms1, atoms2):
    hydrogen_mapping = {}

    for heavy_atom1, heavy_atom2 in zip(heavy_atoms1, heavy_atoms2):
        hydrogen_atoms1 = identify_hydrogen_atoms(heavy_atom1, atoms1)
        hydrogen_atoms2 = identify_hydrogen_atoms(heavy_atom2, atoms2)
        # Assuming hydrogen_atoms1 and hydrogen_atoms2 are lists of Atom objects
        for h1, h2 in zip(hydrogen_atoms1, hydrogen_atoms2):
            hydrogen_mapping[h1.name] = h2.name
    
    return hydrogen_mapping


# Example usage:

# pdb_file1 = 'CLA_602_frank.pdb'
# pdb_file2 = 'CLA_602_mcce.pdb'

# aligned_atoms1, aligned_atoms2 = align_molecules(pdb_file1, pdb_file2)
# print("Molecules aligned successfully.")

# atoms1 = aligned_atoms1
# atoms2 = aligned_atoms2

# if the structure is already aligned, then you can load them directly:
# # atoms1 = load_pdb(pdb_file1)
# # atoms2 = load_pdb(pdb_file2)

# # Identify heavy atoms (non-hydrogen atoms)
# heavy_atoms1 = identify_heavy_atoms(atoms1)
# heavy_atoms2 = identify_heavy_atoms(atoms2)

# # Map heavy atoms based on name and distance
# mapping = map_atoms(heavy_atoms1, heavy_atoms2)
# print("Mapping of heavy atoms based on name and distance:")
# print(mapping)

# # Map hydrogen atoms based on proximity to heavy atoms
# hydrogen_mapping = map_hydrogen_atoms(heavy_atoms1, heavy_atoms2, atoms1, atoms2)
# print("\nMapping of hydrogen atoms based on distance:")
# print(hydrogen_mapping)

# # Combine heavy atom and hydrogen atom mappings into a single dictionary
# combined_mapping = {**mapping, **hydrogen_mapping}
# print("\nCombined mapping of heavy and hydrogen atoms:")
# print(combined_mapping)