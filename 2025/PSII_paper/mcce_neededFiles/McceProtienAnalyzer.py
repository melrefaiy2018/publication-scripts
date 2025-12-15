import matplotlib.pyplot as plt
import seaborn as sns

class McceProteinAnalyzer:
    def __init__(self, pdb_filename):
        self.pdb_filename = pdb_filename
        self.dict_data = {}
        self.list_residue = []
        self.protonation_dict = {}

    def parse_pdb_file(self, list_exclude_residues=None):
        if list_exclude_residues is None:
            list_exclude_residues = ['CLA', 'CHL', 'LHG', 'DGD', 'BCR', 'LHG']
        else:
            list_exclude_residues = list_exclude_residues
            
        with open(self.pdb_filename, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    columns = line.split()
                    residue_name = columns[3]
                    residue_id = columns[4]
                    chain_id = columns[4]
                    charge = columns[-1]
                    residue_identifier = f"{residue_name}_{chain_id}"
                    if residue_identifier.split('_')[0] in list_exclude_residues:
                        continue
                    if charge != 'BK':
                        self.list_residue.append(residue_identifier)
                        self.dict_data[residue_identifier] = charge

    def set_protonation_dict(self, protonation_dict):
        self.protonation_dict = protonation_dict

    def plot_residue_distribution(self, residue_name):
        protonation_states = [charge for residue, charge in self.dict_data.items() if residue.startswith(residue_name)]
        if not protonation_states:
            print(f"No residues found for {residue_name}")
            return
        plt.figure(figsize=(12, 8))
        sns.histplot(protonation_states, bins='auto', kde=False, edgecolor='black', label=residue_name)
        plt.title(f'Protonation Distribution of {residue_name} Residues')
        plt.xlabel('Protonation State')
        plt.ylabel('Frequency')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    def plot_multiple_residue_distributions(self, residue_names):
        plt.figure(figsize=(12, 8))
        for residue_name in residue_names:
            protonation_states = [charge for residue, charge in self.dict_data.items() if residue.startswith(residue_name)]
            if not protonation_states:
                print(f"No residues found for {residue_name}")
                continue
            sns.histplot(protonation_states, bins='auto', kde=False, edgecolor='black', label=residue_name)
        
        plt.title('Protonation Distribution of Selected Residues')
        plt.xlabel('Protonation State')
        plt.ylabel('Frequency')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
        
    # instead of just printin
    # def print_residue_info(self, target_residue):
    #     for resi in self.dict_data.keys():
    #         if target_residue in resi:
    #             print(resi, self.dict_data[resi])

    def get_residue_info(self, target_residue):
        dict_info = {}
        for resi in self.dict_data.keys():
            if target_residue in resi:
                # print(resi, self.dict_data[resi])
                dict_info[resi] = self.dict_data[resi]
        print(dict_info)
        return dict_info

# example usage
# =============
# from pymembrane.util.mcce_util.McceProtienAnalyzer import McceProteinAnalyzer
# # Load the McceProtienAnalyzer object
# mcce_analyzer = McceProteinAnalyzer(f'{md_scratch}{pdb_name}')
# mcce_analyzer.parse_pdb_file()
# mcce_analyzer.plot_residue_distribution('ARG')
# mcce_analyzer.plot_residue_distribution('LYS')
# mcce_analyzer.plot_residue_distribution('ASP')
# mcce_analyzer.plot_residue_distribution('GLU')
# mcce_analyzer.plot_residue_distribution('HIS')
# mcce_analyzer.print_residue_info('LYS')


# # Plot the protonation distribution of LYS residues
# # =================================================
# plt.figure(figsize=(10, 6))
# for residue, charge in mcce_analyzer.dict_data.items():
#     if 'LYS' in residue:
#         plt.bar(residue, charge, color='skyblue', edgecolor='black')
# plt.xlabel('Residue')
# plt.ylabel('Protonation State')
# plt.title('Protonation Distribution of LYS Residues')
# plt.xticks(rotation=45, ha='right')
# plt.grid(axis='y')
# plt.tight_layout()
# plt.show()