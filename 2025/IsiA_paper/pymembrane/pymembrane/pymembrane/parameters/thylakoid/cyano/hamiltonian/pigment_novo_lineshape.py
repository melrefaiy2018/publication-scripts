import numpy as np
from pathlib import Path
from pymembrane.exciton.lineshape import lineshape_drudelorentz, reorganization_drudelorentz,lineshape_brownian_oscillator ,reorganization_brownian_oscillator
from pymembrane.exciton.lineshape import PigmentLineshape_Kubo as Kubo

dir_scratch_wvib = str(Path(__file__).parent.resolve().joinpath('scratch_nova'))
# dir_scratch = str(Path(__file__).parent.resolve().joinpath('scratch_renger_novib'))

w_axis = np.arange(0.01, 2500, 0.01)

cla_brownosc_param = np.loadtxt(
    str(Path(__file__).parent.resolve().joinpath('scratch_nova').joinpath('chlA_BrownianOsc_lhcii.txt')), delimiter=",", dtype=np.float64)

novo_ct = [[lineshape_drudelorentz,     # Overdamped term (Drude-Lorentz)
              reorganization_drudelorentz,
              {'w_axis': w_axis,
               'e_lambda': 200,
               'gamma': 40}],
             [lineshape_brownian_oscillator, # Underdamped term (Brownian Oscillator)
              reorganization_brownian_oscillator,
              {'w_axis': w_axis,
               'e_lambda': cla_brownosc_param[:, 1] * (0.84/0.7004999999999999),
               'gamma': cla_brownosc_param[:, 2],
               'omega': cla_brownosc_param[:, 0]}
              ]]
# result = np.sum(cla_brownosc_param[:, 1] / cla_brownosc_param[:, 0]) = 0.7004999999999999

kubo_novo_cla_ct = Kubo(novo_ct,
                      scratch_dir=dir_scratch_wvib,
                      )
