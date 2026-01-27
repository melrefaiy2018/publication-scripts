import os
import warnings
from Bio import BiopythonWarning, PDB
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBIO, PDBParser
# from pymembrane.util.transformation_matrix import *

warnings.simplefilter('ignore', BiopythonWarning)
plt.rcdefaults()

class ExtendedPDBWriter:
    def __init__(self, filename):
        self.filename = filename

    # add Str methd to this class:
    def __str__(self):
        return f"ExtendedPDBWriter(filename={self.filename})"
    
    # add repr method to this class:
    def __repr__(self):
        return f"ExtendedPDBWriter(filename={self.filename})"
    
    # add eq method to this class:
    def __eq__(self, other):
        return self.filename == other.filename
    

    def write_pdb_line(self, atom_record, **atom_info):
            """
            Writes a line in PDB format to the PDB file.

            Args:
                atom_record (str): The atom record type ("ATOM" or "HETATM").
                **atom_info: Keyword arguments containing atom information.
                    - residue_name (str): The name of the residue.
                    - atom_index (int): The index of the atom.
                    - atom_name (str): The name of the atom.
                    - chain_id (str): The chain identifier.
                    - residue_index (int): The index of the residue.
                    - x (float): The x-coordinate of the atom.
                    - y (float): The y-coordinate of the atom.
                    - z (float): The z-coordinate of the atom.
                    - element_name (str): The name of the element.
                    - charge (int): The charge of the atom.

            Returns:
                None
            """
            # print(f'{atom_info=}')
            if atom_info['residue_name'] in ("CLA", "CLB", "CHL", "BCR", "SQD", 'LHG', 'LMG', 'LMU', 'PQN','DGD'):
                atom_record = "HETATM"
            else:
                atom_record = "ATOM"

            pdb_line = self._format_atom_record(atom_record)
            pdb_line += self._format_atom_index(atom_info['atom_index'])
            pdb_line += self._format_atom_name(atom_info['atom_name'])
            pdb_line += self._format_residue_name(atom_info['residue_name'])
            pdb_line += self._format_chain_id(atom_info['chain_id'])
            pdb_line += self._format_residue_index(atom_info['residue_index'])
            pdb_line += self._format_xyz(atom_info['x'])
            pdb_line += self._format_xyz(atom_info['y'])
            pdb_line += ' '
            pdb_line += self._format_xyz(atom_info['z'])
            pdb_line += self._format_element_name(atom_info['element_name'])
            pdb_line += self._format_charge(atom_info['residue_name'], atom_info['charge'])
            pdb_line += "\n"

            with open(self.filename, "a") as pdb_file:
                pdb_file.write(pdb_line)

    def _format_atom_record(self, atom_record):
        return f"{atom_record:<6s}"

    def _format_atom_index(self, atom_index):
        return f"{atom_index:>5d} "

    def _format_atom_name(self, atom_name):
        return f"{self._construct_atomtype_string(atom_name)} "

    def _format_residue_name(self, residue_name):
        return f"{residue_name:>3s} "

    def _format_chain_id(self, chain_id):
        return f"{chain_id:>1s}"

    def _format_residue_index(self, residue_index):
        return f"{residue_index:>4s}    "

    def _format_xyz(self, my_numb):
        my_numb = float(my_numb)
        return f'{my_numb:>8.3f}'

    def _format_element_name(self, element_name):
        return " " * 21 + f"{element_name:>2s}  "

    def _format_charge(self, residue_name, charge):
        # print(f'{residue_name=}, {charge=}')
        if residue_name in ("CLA", "CLB", 'CHL'):
            return "  None"
        else:
            return f"{charge:7.4f}"

    def _construct_atomtype_string(self, my_atomtype):
        # 13-16 columns (4 columns)
        if len(my_atomtype) > 4:
            raise ValueError(f"atom type can't be me more than 5 characters")
        elif len(my_atomtype) == 4:
            return f'{my_atomtype:<4s}'

        elif len(my_atomtype) == 3:
            return f'{my_atomtype:>4s}'
        else:
            return f'{my_atomtype:^4s}'

    def write_pdb_file(self, protein_atomic):
        for chain in protein_atomic._pdb.get_chains():
            for residue in chain.get_residues():
                for atom in residue.get_atoms():
                    if atom.get_name() == 'MG':
                        atom_element = 'Mg'
                    else:
                        atom_element = atom.element

                    atom_info = {
                        "atom_record": "ATOM",
                        "atom_index": atom.get_serial_number(),
                        "atom_name": atom.get_name(),
                        "residue_name": residue.resname,
                        "chain_id": chain.id,
                        "residue_index": str(residue.id[1]),
                        "x": float(atom.coord[0]),
                        "y": float(atom.coord[1]),
                        "z": float(atom.coord[2]),
                        "element_name": atom_element,
                        "charge": atom.get_charge() if atom.get_charge() is not None else 0,
                    }
                    self.write_pdb_line(**atom_info)




# # Example run:
# # ============
# path = ''
# mcce_pdb = 'extended_CP29_pH8_chain_r02.pdb' # mcce extended pdb that has charges
# from pymembrane.structure.atomic_protein import ElectrostaticProteinAtomic

# pdb_atomic = PDB.PDBParser().get_structure('LHCII', mcce_pdb)
# lhcii_atomic = ElectrostaticProteinAtomic(pdb_atomic, path_extended_pdb=mcce_pdb, name='LHCII')


# # you can do some transformation here on the protein_atomic_1
# # Now, save the transformed protein_atomic_1 to a new pdb file called mcce_pdb_oriented.pdb:
# output_pdb = 'extended_CP29_pH8_chain_r02.pdb'
# writer = ExtendedPDBWriter(output_pdb)
# writer.write_pdb_file(lhcii_atomic)


# pdb_atomic_out = PDB.PDBParser().get_structure('LHCII', output_pdb)
# lhcii_atomic_out = ElectrostaticProteinAtomic(pdb_atomic_out, path_extended_pdb=output_pdb, name='LHCII')

# for chain in lhcii_atomic._pdb.get_chains():
#     for residue in chain.get_residues():
#         for atom in residue.get_atoms():
#             print(atom.get_name(), atom.get_coord())

