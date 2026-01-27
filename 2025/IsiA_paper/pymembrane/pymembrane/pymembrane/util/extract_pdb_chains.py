from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser
import os


class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:
            return 0


def extract_chain_from_pdb(PDB_name, PDB_path, list_chains):
    structure = PDBParser(QUIET=True, PERMISSIVE=1).get_structure(PDB_name, PDB_path)
    print(structure)
    for chain in list_chains:
        pdb_chain_file = str(structure.id) + '_chain_{}.pdb'.format(chain)
        io_w_no_h = PDBIO()
        io_w_no_h.set_structure(structure)
        io_w_no_h.save('{}'.format(pdb_chain_file), ChainSelect(chain))
        print(os.getcwd())



# example usage:
# =============================================================================
# list_chains = ['A','D']
# PDB_path = f"{os.getcwd()}/3ARC.pdb" # path of the input pdb file that you want to extract the chains form.
# PDB_name = '3ARC' # prefex of the output file
# extract_chain_from_pdb(PDB_name, PDB_path, list_chains)


