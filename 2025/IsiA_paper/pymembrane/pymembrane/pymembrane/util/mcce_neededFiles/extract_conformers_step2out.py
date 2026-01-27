# this script will look for the occupancy of the conformer id you provide by looking
# into fort.38 file and then extract their coordinates as a pdb file by looking for
# their ids in step2_out.pdb file.

# Mohamed Elrefaiy
# Aug.19, 2022


working_dir = '/Users/mohamed/Documents/Research/Projects/WSCP/2024_reppert/new_run/run_with_water/Q57D/'

step2_out = working_dir + 'step2_out.pdb'
fort38 = working_dir + 'fort.38'

# searching_name = 'LHG01R2630'
searching_name = 'ASP02A0057'
searching_name = 'ASP01A0057'
searching_name = 'ASP-1D0057'


lst_conf_ids_occupancy = []

for line in open(fort38).readlines()[1:]:
    # print(line[0:10])
    if str(searching_name).strip() == line[0:10].strip():
        # line[15:20] is the occupancy:
        if float(line[15:20]) > float(0.000):
            dict_conf_ids_occupancy = {'conf_id': line[0:14], 'occupancy': line[15:20]}
            lst_conf_ids_occupancy.append(dict_conf_ids_occupancy)

print('lst_conf_ids_occupancy = ', lst_conf_ids_occupancy)
# build directory called ddd:
import os
ddd = working_dir + 'extract_confs'
if not os.path.exists(ddd):
    os.makedirs(ddd)

for index in lst_conf_ids_occupancy:
    conf_id = index['conf_id']  # LHG01R2630_275
    with open(f'{working_dir}extract_confs/{conf_id}_{index["occupancy"]}.pdb',
              'w') as conf_pdb:
        for line in open(step2_out).readlines():
            res_name = line[17:20]  # LHG
            conf_name = line[21:26]  # R2630
            conf_number = line[27:30]  # 001
            conf_charge = line[80:82]  # 01 or BK

            if conf_name == conf_id[5:10]:  # conf_id[5:10] = R0099
                if conf_charge == 'BK':
                    extended_line = f'{line[0:57]}'"                   "f'' \
                                    f'{line[57:90]}\n'
                    conf_pdb.write(extended_line)
                else:
                    if conf_id[11:15] == conf_number:  # conf_id[11:15] = 002
                        extended_line = f'{line[0:57]}'"                   "f'' \
                                        f'{line[57:90]}\n'
                        conf_pdb.write(extended_line)
