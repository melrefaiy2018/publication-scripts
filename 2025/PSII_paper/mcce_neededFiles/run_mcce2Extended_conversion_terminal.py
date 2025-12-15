#!/usr/bin/env python
# to install this script in the enviroment and run it from the terminal, use the following commands:
# conda info --envs > this will give you the path of all installed conda enviroments, take the path of the enviroment you want to install the script in.
# mkdir -p <eniroment path>/bin

# Create the Symbolic Link:
# ln -s $(pwd)/run_mcce2Extended_conversion_terminal.py <eniroment path>/bin/

# Now, in the terminal, you can run the script from any directory using the following command:
# conda activate pymembrane_alpha
# run_mcce2Extended_conversion_terminal.py -i most_occ_pH8.pdb -o extended_most_occ_pH8.pdb -d MEM

import argparse
from pymembrane.util.mcce_neededFiles.mcce2extendedPDB_converter import Mcce2ExtendedConverter

def main(args):
    pdb_file_input = args.input_pdb
    pdb_file_output = args.output_pdb
    list_ligand_delete = args.list_ligand_2_delete.split(',')  # assuming comma-separated list of ligands
    converter = Mcce2ExtendedConverter(mcce_pdb_path=pdb_file_input,
                                        list_ligand_delete=list_ligand_delete,
                                        custom_charge=args.custom_charge)
    converter.convert_mcce_to_pymembrane(pdb_output_name=pdb_file_output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MCCE to Extended PDB Converter')
    parser.add_argument('-i', '--input_pdb', help='Input PDB file', type=str, required=True)
    parser.add_argument('-o', '--output_pdb', help='Output PDB file', type=str, required=True)
    parser.add_argument('-d', '--list_ligand_2_delete', help='Comma-separated list of ligands to be deleted', type=str, required=True)
    parser.add_argument('-c', '--custom_charge', help='Use custom charges', action='store_true')
    args = parser.parse_args()
    main(args)
    