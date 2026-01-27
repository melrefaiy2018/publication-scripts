from testing.test_atomic_protein import FakeAtom, FakeResidue

from pymembrane.structure.atomic_pigment import *
from pymembrane.util.physical_constants import precision

# Chl
atom0 = FakeAtom((0, 0, 0), mass=14, element='N', id=0)  # NA
atom1 = FakeAtom((0, 2, 0), mass=14, element='N', id=1)  # ND
atom2 = FakeAtom((2, 0, 0), mass=14, element='N', id=2)  # NB
atom3 = FakeAtom((2, 2, 0), mass=14, element='N', id=3)  # NC
atom4 = FakeAtom((1, 1, 0), mass=24, element='MG', id=4)
dict_atoms_1 = {'NA': atom0, 'NB': atom2, 'NC': atom3, 'ND': atom1, 'MG': atom4}

# Cla
atom5 = FakeAtom((5, 5, 5), mass=14, element='N', id=5)  # NA
atom6 = FakeAtom((5, 7, 5), mass=14, element='N', id=6)  # ND
atom7 = FakeAtom((7, 5, 5), mass=14, element='N', id=7)  # NB
atom8 = FakeAtom((7, 7, 5), mass=14, element='N', id=8)  # NC
atom9 = FakeAtom((6, 6, 5), mass=24, element='MG', id=9)
dict_atoms_2 = {'NA': atom5, 'NB': atom7, 'NC': atom8, 'ND': atom6, 'MG': atom9}

# build fake residue
residue_1 = FakeResidue(dict_atoms_1)
residue_2 = FakeResidue(dict_atoms_2)

# add data to the fake residues:
chl = ChlorophyllAtomic(residue_1)
chl_tresp_dict = {'NA': -0.026, 'NB': 0.086, 'NC': 0.041, 'ND': -0.078, 'MG': 0.0040}
chl.dict_data['tresp_atoms'] = list(chl_tresp_dict.keys())
chl.dict_data['tresp_pc'] = list(chl_tresp_dict.values())
chl.dict_data['dipole_mag'] = 4.6
chl.dict_data['vacuum_mag'] = 1
chl.dict_data['atom_pos'] = list([np.array(chl.residue.dict_atoms[k].coord) for k in chl.residue.dict_atoms])

cla = ChlorophyllAtomic(residue_2)
cla_tresp_dict = {'NA': 0.03, 'NB': -0.062, 'ND': 0.12, 'NC': -0.012, 'MG': -0.022}
cla.dict_data['tresp_atoms'] = list(cla_tresp_dict.keys())
cla.dict_data['tresp_pc'] = list(cla_tresp_dict.values())
cla.dict_data['dipole_mag'] = 4.4
cla.dict_data['vacuum_mag'] = 1
cla.dict_data['atom_pos'] = list([np.array(cla.residue.dict_atoms[k].coord) for k in cla.residue.dict_atoms])


def test_scale_charges():
    # CLA case:
    assert np.isclose(scale_charges(cla.dict_data['tresp_pc'], cla.dict_data['atom_pos'],
                                    cla.dict_data['vacuum_mag']), 0.648607, atol=precision)
    # CHL case:
    assert np.isclose(scale_charges(chl.dict_data['tresp_pc'], chl.dict_data['atom_pos'], chl.dict_data['vacuum_mag']),
                      0.2673275, atol=precision)

def test_calculate_tresp_coupling():
    # Testing using defalt values for the tresp_CC = 1.1615E5 (units of Angstrom*cm^-1)/e^2
    assert np.isclose(calculate_tresp_coupling(chl, cla, f_val=1, dielectric=1), 128.41959475, atol=precision)

    # testing using non-dafalt values:
    assert np.isclose(calculate_tresp_coupling(chl, cla, f_val=10, dielectric=1, coupling_const=1), 0.011056357705, atol=precision)
    assert np.isclose(calculate_tresp_coupling(chl, cla, f_val=1, dielectric=20, coupling_const=1), 5.52817e-05, atol=precision)
    assert np.isclose(calculate_tresp_coupling(chl, cla, f_val=1, dielectric=1, coupling_const=30), 0.0331691, atol=precision)

def test_calculate_dipole_coupling():
    # Testing using defalt values for the dipole_CC = (1 / (4 * np.pi * 1.5812E-5)) ,units of (cm^-1 *Angstrom^3)/D^2:
    assert np.isclose(calculate_dipole_coupling(chl, cla, dielectric=1), 156.82740033360128, atol=precision)

    # testing using non-dafalt values:
    assert np.isclose(calculate_dipole_coupling(chl, cla, dielectric=1, coupling_const=1), 0.0311615, atol=precision)
    assert np.isclose(calculate_dipole_coupling(chl, cla, dielectric=1, coupling_const=20), 0.623230, atol=precision)
    assert np.isclose(calculate_dipole_coupling(chl, cla, dielectric=200, coupling_const=1), 0.00015580, atol=precision)

