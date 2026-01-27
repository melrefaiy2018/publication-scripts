import numpy as np
from pathlib import Path
from pymembrane.exciton.lineshape import lineshape_renger, reorganization_renger
from pymembrane.exciton.lineshape import PigmentLineshape_Kubo as Kubo

dir_scratch_wvib = str(Path(__file__).parent.resolve().joinpath('scratch_renger'))
dir_scratch = str(Path(__file__).parent.resolve().joinpath('scratch_renger_novib'))

list_w_vib = [100, 175, 250, 300, 375, 500, 600, 725, 800, 875]
list_s_vib = [0.2, 0.1, 0.06, 0.04, 0.06, 0.04, 0.015, 0.04, 0.02, 0.02]

# This list for the IsiA monomer calculation (CT):
list_s_vib_reduced = [0.2, 0.1, 0.06, 0.04, 0.06, 0.04, 0, 0, 0, 0]

lineshape_parameters = [
    [lineshape_renger,
     reorganization_renger,
     {'w_axis': np.arange(0.01, 2000, 0.01),
      'S0': 0.5,
      's1': 0.8,
      's2': 0.5,
      'w1': 0.069 * 8.0656,  # Units: cm^-1
      'w2': 0.24 * 8.0656  # Units: cm^-1
      }
     ]]
# lineshape for charge transfer state (CT):
# ------------------------------------------
CT_lineshape_parameters = [
    [lineshape_renger,
     reorganization_renger,
     {'w_axis': np.arange(0.01, 2000, 0.01),
      'S0': 0.5,
      's1': 0.8,
      's2': 0.5,
      'w1': 0.069 * 8.0656,  # Units: cm^-1
      'w2': 0.24 * 8.0656  # Units: cm^-1
      }
     ]]

kubo_renger_cla = Kubo(lineshape_parameters,
                       explicit_vib=[list_s_vib, list_w_vib],
                       scratch_dir=dir_scratch_wvib,
                       )

kubo_renger_rc = Kubo(lineshape_parameters,
                      scratch_dir=dir_scratch,
                      )

kubo_renger_cla_ct = Kubo(CT_lineshape_parameters,
                      scratch_dir=dir_scratch,
                      )
kubo_renger_cla_reduce = Kubo(lineshape_parameters,
                       explicit_vib=[list_s_vib_reduced, list_w_vib],
                       scratch_dir=dir_scratch_wvib,
                       )