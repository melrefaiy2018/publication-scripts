import os
import numpy as np
from pymembrane.util.mcce_neededFiles.psii_cofactor_charges_frank import dict_LUT_frank, dict_NEX_frank, dict_XAT_frank
from collections import defaultdict, Counter

list_frank_charges = [dict_LUT_frank, 
                      dict_NEX_frank, 
                    #   dict_XAT_frank
                      ]

class Mcce2ExtendedConverter:
    """
    A class for converting MCCE PDB files to extended PDB format.
    """
    standard_amino_acids = {
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL','CYD'
            }

    def __init__(self, mcce_pdb_path, list_ligand_delete=None, custom_charge=False):
        """
        Initialize the class with the MC
        Parameters:
        ----------
        1.) mcce_pdb_path (str): The path to the MCCE PDB file.
        2.) path_saving_dir (str): The directory to save the converted files. If None, the current working directory is used.

        Returns:
        --------
        1.) None

        """
        self.mcce_pdb_path = mcce_pdb_path
        self.list_ligand_delete = list_ligand_delete
        self.path_saving_dir = os.getcwd()
        self.custom_charge = custom_charge  # Add custom_charge flag to the class

        
    def _delete_intermediate_files(self, file_list):
        """
        Delete the intermediate files from the given file list.

        Parameters:
            file_list (list): List of file paths to be deleted.

        Returns:
            None
        """
        for file_path in file_list:
            if os.path.exists(file_path):
                os.remove(file_path)

    def convert_mcce_to_pymembrane(self, pdb_output_name='extended_most_occ_mcce.pdb'):
        # 1. Convert mcce PDB format to extended PDB format:
        print(f'Converting MCCE PDB to extended PDB...\n')
        self.convert_to_extended_pdb_format()

        # 2. Rename the N and C terminal caps that mcce assigns to the terminal residues:
        print(f'Renaming capping residues...\n')
        self.rename_capping_residues(self.extended_pdb_1)
        final_pdb_name = os.path.join(self.path_saving_dir, pdb_output_name)

        # 3. Delete ligands:
        if self.list_ligand_delete is not None:
            print(f'Deleting ligands...\n')
            self.delete_ligand(self.extended_pdb_2, final_pdb_name, self.list_ligand_delete)

        else:
            os.rename(self.extended_pdb_2, final_pdb_name)

        # Delete intermediate files
        file_list_delete = [self.extended_pdb_1, self.extended_pdb_2]
        print(f'Deleting intermediate files...\n')
        self._delete_intermediate_files(file_list_delete)

        print(f'Conversion is done!\n')
        print(f'PDB saved in dir : {self.path_saving_dir}\n')


    def convert_to_extended_pdb_format(self):  # fix the input file, by calling it from the init method.
        """
        Converts an MCCE PDB file format to an extended PDB file format.

        Parameters:
            mcce_pdb_path (str): Path to the MCCE PDB file.
            extended_pdb_path_0 (str): Path to the extended PDB file.

        Returns:
            None
        """
        self.extended_pdb_1 = os.path.join(self.path_saving_dir, "extended_pdb1.pdb")
        with open(self.mcce_pdb_path, 'r') as mcce_pdb, open(self.extended_pdb_1, 'w') as extended_pdb:
            for mcce_line in mcce_pdb:
                if mcce_line.startswith('ATOM') or mcce_line.startswith('HETATM'):
                    extended_line = self.__construct_extended_pdb_line(mcce_line)
                    extended_pdb.write(extended_line)

    
    def __construct_extended_pdb_line(self, mcce_line):
        """
        Construct the extended PDB line from the MCCE PDB line.
        """
        atom_info = self.parse_pdb_line_string(mcce_line)
        residue_name = atom_info['residue_name']
        atom_type = atom_info['atom_type']
        chain_id = atom_info['chain_id']
        element_symbol = atom_info['element_symbol']
        mcce_conf_number = atom_info['mcce_conf_number']
        mcce_charges = atom_info['mcce_charges']

        # Dictionary mapping residue names to their respective charge dictionaries
        charge_dicts = {
            'LUT': dict_LUT_frank,
            'NEX': dict_NEX_frank,
        }

        # Decide on the charge
        if self.custom_charge and residue_name in charge_dicts and atom_type in charge_dicts[residue_name]:
            charge = self.__construct_charge_string(charge_dicts[residue_name][atom_type])
        elif residue_name in ['CLA', 'CLB', 'CHL', 'CT1', 'CT2', 'CT3', 'CT4']:
            charge = 'None'  # Special case for ligands
        else:
            charge = mcce_charges  # Use MCCE charges if custom_charge is False or no dictionary entry exists

        # Determine the record type
        record_name = 'ATOM' if residue_name in self.standard_amino_acids else 'HETATM'

        # Start constructing the base PDB line
        base_line = (
            self.__construct_record_name_string(record_name) +
            self.__construct_atom_serial_number_string(atom_info['atom_index']) + " " +
            self.__construct_atomtype_string(atom_type) +
            self.__construct_altloc_string() +
            self.__construct_residue_name_string(residue_name) + " "+
            self.__construct_chain_id(chain_id) +
            self.__construct_residue_seq_string(atom_info['residue_index']) +
            self.__construct_icode_string() + "   " +
            self.__construct_xyz_string(atom_info['x_coord']) +
            self.__construct_xyz_string(atom_info['y_coord']) +
            self.__construct_xyz_string(atom_info['z_coord']) + " "
        )
        # Construct the complete extended line
        extended_line = (
            base_line + ' ' * 21 +
            f"{self.__construct_element_symbol_string(element_symbol)}  " +
            f"{charge:>6} {mcce_conf_number:>6} \n"
        )

        return extended_line

    # @staticmethod
    def parse_pdb_line_string(self, pdb_line):
        """
        Parse the PDB line string and return the atom information in a dictionary format.
        """
        atom_info = {
            'record_atom': pdb_line[0:6].strip(),
            'atom_index': pdb_line[6:11].strip(),
            'atom_type': pdb_line[12:16].strip(),
            'residue_name': pdb_line[17:20].strip(),
            'chain_id': pdb_line[21:22].strip(),
            'residue_index': pdb_line[22:26].strip(),
            'x_coord': pdb_line[30:38].strip(),
            'y_coord': pdb_line[38:46].strip(),
            'z_coord': pdb_line[46:54].strip(),
            'element_symbol': pdb_line[76:78].strip(),
            'mcce_charges': pdb_line[68:76].strip(),  
            'mcce_conf_number': pdb_line[80:82].strip()
        }

        # if not atom_info['element_symbol'] and atom_info['atom_type'] == 'MG':
        if atom_info['atom_type'] == 'MG':
            atom_info['element_symbol'] = 'Mg'
        elif not atom_info['element_symbol']:
            atom_info['element_symbol'] = atom_info['atom_type'][0]

        return atom_info


    def rename_capping_residues(self, extended_pdb_path_0):
            """
            Rename the capping residues in the extended PDB file, except for HOH.
            """
            self.extended_pdb_2 = os.path.join(self.path_saving_dir, 'extended_pdb2.pdb')

            # Step 1: Read the PDB file and build a mapping of residue indices to their corresponding residue names.
            residue_map = defaultdict(list)
            with open(extended_pdb_path_0, 'r') as extended_pdb:
                for line in extended_pdb:
                    residue_index = line[21:26].strip()
                    residue_name = line[17:20].strip()
                    residue_map[residue_index].append(residue_name)

            # Step 2: Determine the most frequent residue name for each index, prioritizing uppercase names.
            unified_residue_map = {}
            for index, names in residue_map.items():
                # Skip HOH residues in the mapping process
                if 'HOH' in names:
                    continue

                name_counts = Counter(names)
                uppercase_names = {name: count for name, count in name_counts.items() if name.isupper()}
                lowercase_names = {name: count for name, count in name_counts.items() if not name.isupper()}

                if uppercase_names:
                    most_common_name = max(uppercase_names, key=uppercase_names.get)
                else:
                    most_common_name = max(lowercase_names, key=lowercase_names.get)

                unified_residue_map[index] = most_common_name

            # Step 3: Rewrite the PDB file with unified residue names based on the mapping, skipping HOH.
            with open(extended_pdb_path_0, 'r') as extended_pdb, open(self.extended_pdb_2, 'w') as complex_pdb:
                for line in extended_pdb:
                    residue_index = line[21:26].strip()
                    residue_name = line[17:20].strip()
                    
                    # Skip HOH lines by checking for both HETATM and HOH
                    if line.startswith("HETATM") and residue_name == "HOH":
                        complex_pdb.write(line)
                        continue

                    # Use the unified name if it exists in the mapping
                    unified_resname = unified_residue_map.get(residue_index, residue_name)
                    complex_line = line[:17] + unified_resname.ljust(3) + line[20:]
                    complex_pdb.write(complex_line)


    def delete_ligand(self, extended_pdb_path_0, extended_pdb_path_1, list_ligand_delete):
        """
        Delete the ligands from the extended PDB file.
        
        Parameters:
        -----------
        1.) extended_pdb_path_0 (str): The path to the extended PDB file.
        2.) extended_pdb_path_1 (str): The path to the complex PDB file.
        3.) list_ligand_delete (list): List of ligands to delete.
        
        Returns:
        --------
        1.) None

        """
        with open(extended_pdb_path_0, 'r') as extended_pdb, open(extended_pdb_path_1, 'w') as complex_pdb:
            for extended_line in extended_pdb:
                residue_name = extended_line[17:20].strip()

                if residue_name in list_ligand_delete:
                    pass

                else:
                    complex_pdb.write(extended_line)

    @staticmethod
    def __construct_record_name_string(atom_record):
        """
        Construct the record name string.
        
        Parameters:
        -----------
        atom_record (str): The atom record name.
        
        Returns:
        --------
        str: The formatted record name string.
        """
        return f'{atom_record:<6}'

    @staticmethod
    def __construct_atom_serial_number_string(atom_index):
        """
        Construct the atom serial number string.
        
        Parameters:
        -----------
        atom_index (str): The atom index.

        Returns:
        --------
        str: The formatted atom index string.
        """
        return f'{str(atom_index):>5}'

    @staticmethod
    def __construct_atomtype_string(my_atomtype):
        """
        Construct the atom type string.

        Parameters:
        -----------
        my_atomtype (str): The atom type.

        Returns:
        --------
        str: The formatted atom type string.
        """
        if len(my_atomtype) > 4:
            raise ValueError("Atom type can't be more than 4 characters")
        return f'{my_atomtype:>4}'

    @staticmethod
    def __construct_altloc_string():
        """
        Construct the alternate location indicator string.

        Returns:
        --------
        str: The formatted alternate location indicator string.
        """
        return ' '

    @staticmethod
    def __construct_residue_name_string(residue_name):
        """
        Construct the residue name string.
        
        Parameters:
        -----------
        residue_name (str): The residue name.

        Returns:
        --------
        str: The formatted residue name string.
        """
        return f'{residue_name:<3}'

    @staticmethod
    def __construct_residue_seq_string(residue_index):
        """
        Construct the residue sequence number string.
        
        Parameters:
        -----------
        residue_index (str): The residue sequence number.

        Returns:
        --------
        str: The formatted residue sequence number string.
        """
        return f'{residue_index:>4}'

    @staticmethod
    def __construct_icode_string():
        """
        Construct the insertion code string.

        Returns:
        --------
        str: The formatted insertion code string.
        """
        return ' '

    @staticmethod
    def __construct_xyz_string(my_numb):
        """
        Construct the XYZ coordinate string.

        Parameters:
        -----------
        my_numb (str): The number to be converted to a string.

        Returns:
        --------
        str: The formatted string.
        """
        my_numb = float(my_numb)
        return f'{my_numb:>8.3f}'

    @staticmethod
    def __construct_element_symbol_string(element_name):
        """
        Construct the element symbol string.

        Parameters:
        -----------
        element_name (str): The element name.

        Returns:
        --------
        str: The element name string.
        """
        return f'{element_name:<2}'

    @staticmethod
    def __construct_charge_string(charge):
        """
        Construct the charge string.

        Parameters:
        -----------
        charge (str): The charge.

        Returns:
        --------
        str: The formatted charge string.
        """
        return f'{float(charge):6.3f}'
    
    @staticmethod
    def __construct_chain_id(chain_id):
        """
        Construct the element symbol string.

        Parameters:
        -----------
        element_name (str): The element name.

        Returns:
        --------
        str: The element name string.
        """
        return f'{chain_id:<1}'