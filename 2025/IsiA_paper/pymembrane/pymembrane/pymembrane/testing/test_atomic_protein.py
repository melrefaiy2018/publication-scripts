import os

import numpy as np
import pandas as pd
import pytest

from pymembrane.structure.atomic_pigment import ChlorophyllAtomic
from pymembrane.structure.atomic_protein import ProteinAtomic
from pymembrane.util.fake_objects_for_testing import FakePDB, FakeChild, FakeAtom, FakeResidue
from pymembrane.util.physical_constants import precision


@pytest.fixture(scope="function")
def linear_rotors():
    H2_pdb = FakePDB([FakeAtom([-0.5, -0.5, -0.5]), FakeAtom([0.5, 0.5, 0.5])], id=1)
    H2 = ProteinAtomic(pdb=H2_pdb, name="H2")
    H2_com = np.array([0, 0, 0], dtype=float)
    H2_I = np.array([[1, -0.5, -0.5], [-0.5, 1, -0.5], [-0.5, -0.5, 1]], dtype=float)

    HD_pdb = FakePDB([FakeAtom([-2 / 3, -2 / 3, -2 / 3]), FakeAtom([1 - 2 / 3, 1 - 2 / 3, 1 - 2 / 3], mass=2)])
    HD = ProteinAtomic(pdb=HD_pdb, name="HD")
    HD_com = np.array([0, 0, 0], dtype=float)
    HD_I = [[1.33333333, -0.66666667, -0.66666667], [-0.66666667, 1.33333333, -0.66666667],
            [-0.66666667, -0.66666667, 1.33333333]]

    yield {"H2_pdb": H2_pdb, "H2": H2, "H2_com": H2_com, "H2_I": H2_I, "HD_pdb": HD_pdb, "HD": HD, "HD_com": HD_com,
           "HD_I": HD_I}


@pytest.fixture(scope="function")
def spherical_rotors():
    C = FakeAtom((0, 0, 0), mass=12)
    Ha = FakeAtom((np.sqrt(8 / 9), 0, -1 / 3))
    Hb = FakeAtom((-np.sqrt(2 / 9), np.sqrt(2 / 3), -1 / 3))
    Hc = FakeAtom((-np.sqrt(2 / 9), - np.sqrt(2 / 3), -1 / 3))
    Hd = FakeAtom((0, 0, 1))
    CH4_pdb = FakePDB([C, Ha, Hb, Hc, Hd], id=1)
    CH4 = ProteinAtomic(pdb=CH4_pdb, name="CH4")
    CH4_I = [[2.66666667, 0, 0], [0, 2.66666667, 0], [0, 0, 2.66666667]]
    CH4_com = np.array([0, 0, 0], dtype=float)
    yield {"C": C, "Ha": Ha, "Hb": Hb, "Hc": Hc, "Hd": Hd, "CH4_pdb": CH4_pdb, "CH4": CH4, "CH4_I": CH4_I,
           "CH4_com": CH4_com}


@pytest.fixture(scope="function")
def prolate_rotors():
    Cl = FakeAtom((0, 0, 1 - 0.68), mass=35)
    C = FakeAtom((0, 0, -0.68), mass=12)
    Ha = FakeAtom((np.sqrt(8 / 9), 0, -1 / 3 - 0.68))
    Hb = FakeAtom((-np.sqrt(2 / 9), np.sqrt(2 / 3), -1 / 3 - 0.68))
    Hc = FakeAtom((-np.sqrt(2 / 9), -np.sqrt(2 / 3), -1 / 3 - 0.68))
    CH3Cl_pdb = FakePDB([C, Ha, Hb, Hc, Cl], id=1)
    CH3Cl = ProteinAtomic(pdb=CH3Cl_pdb, name="CH3Cl")
    CH3Cl_com = np.array([0, 0, 0])
    CH3Cl_I = [[36.66666666666667, 0, 0], [0, 36.66666666666667, 0], [0, 0, 2.6666666666666665]]

    yield {"Cl": Cl, "C": C, "Ha": Ha, "Hb": Hb, "Hc": Hc, "CH3Cl_pdb": CH3Cl_pdb, "CH3Cl": CH3Cl,
           "CH3Cl_com": CH3Cl_com, "CH3Cl_I": CH3Cl_I}


@pytest.fixture(scope="function")
def asymmetric_rotors():
    Ha = FakeAtom((np.sqrt(8 / 9) - 0.02618914, -0.045360921162651446, -1 / 3 - 0.037037037037037035))
    Hb = FakeAtom((-np.sqrt(2 / 9) - 0.02618914, np.sqrt(2 / 3) - 0.045360921162651446, -1 / 3 - 0.037037037037037035))
    O = FakeAtom((-0.02618914, -0.045360921162651446, - 0.037037037037037035), mass=16)
    H2O_pdb = FakePDB([O, Ha, Hb])
    H2O = ProteinAtomic(pdb=H2O_pdb, name="H2O")
    H2O_I = [[0.82716049, 0.40628352, 0.13967541], [0.40628352, 1.2962963, 0.24192491],
             [0.13967541, 0.24192491, 1.72839506]]

    H2O_mass = Ha.mass + Hb.mass + O.mass

    H2O_com = np.zeros(3)  # this is not the actual com, but it is expected that the pdb is shifted after inizialization

    H2O_shifted = ProteinAtomic(pdb=H2O_pdb, name="H2O", center_atomic=False)
    H2O_shifted_com = (Ha.coord * Ha.mass + Hb.coord * Hb.mass + O.coord * O.mass) / H2O_mass  # THIS is the actual com
    H2O_shifted_I = H2O_I

    yield {"O": O, "H2O_pdb": H2O_pdb, "H2O": H2O, "H2O_I": H2O_I, "H2O_com": H2O_com, "H2O_mass": H2O_mass,
           "H2O_name": 'H2O', "H2O_shifted": H2O_shifted, "H2O_shifted_pdb": H2O_shifted,
           "H2O_shifted_com": H2O_shifted_com, "H2O_shifted_I": H2O_shifted_I}


# TESTING SUITE
# -------------
@pytest.mark.parametrize("rotor_type, atomic_prot, atomic_prot_pdb, atomic_prot_com",
                         [
                             ("linear_rotors", "H2", "H2_pdb", "H2_com"),
                             ("linear_rotors", "HD", "HD_pdb", "HD_com"),
                             ("spherical_rotors", "CH4", "CH4_pdb", "CH4_com"),
                             ("prolate_rotors", "CH3Cl", "CH3Cl_pdb", "CH3Cl_com"),
                             ("asymmetric_rotors", "H2O", "H2O_pdb", "H2O_com"),
                             ("asymmetric_rotors", "H2O_shifted", "H2O_pdb", "H2O_shifted_com")
                         ])
def test_initialization(rotor_type, atomic_prot, atomic_prot_pdb, atomic_prot_com, request):
    fixture = request.getfixturevalue(rotor_type)
    assert fixture[atomic_prot].pdb == fixture[atomic_prot_pdb]
    assert np.allclose(fixture[atomic_prot].center_of_mass, fixture[atomic_prot_com], atol=precision)
    assert fixture[atomic_prot].dict_data == {}
    assert fixture[atomic_prot].dict_pigments == {}


@pytest.mark.parametrize("rotor_type, atomic_protein, same_atomic_protein, different_atomic_protein",
                         [
                             ("linear_rotors", "H2", "H2", "HD")
                         ])
def test_eq(rotor_type, atomic_protein, same_atomic_protein, different_atomic_protein, request):
    fixture = request.getfixturevalue(rotor_type)
    assert fixture[atomic_protein].name == fixture[same_atomic_protein].name != fixture[different_atomic_protein].name


@pytest.mark.parametrize("rotor_type, atomic_protein, I",
                         [("linear_rotors", "H2", "H2_I"),
                          ("linear_rotors", "HD", "HD_I"),
                          ("spherical_rotors", "CH4", "CH4_I"),
                          ("asymmetric_rotors", "H2O", "H2O_I"),
                          ("asymmetric_rotors", "H2O_shifted", "H2O_shifted_I"),
                          ])
def test_inertia_tensor(rotor_type, atomic_protein, I, request):
    fixture = request.getfixturevalue(rotor_type)
    assert np.allclose(fixture[atomic_protein].inertia_tensor, fixture[I], atol=precision)


@pytest.mark.parametrize("rotor_type, atomic_protein, atom_coord",
                         [
                             ("linear_rotors", "H2", [[-0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
                             ("linear_rotors", "HD",
                              [[-0.66666669, -0.66666669, 0.66666665], [0.33333331, 0.33333331, -0.33333335]]),
                             ("spherical_rotors", "CH4",
                              [[0, 0, 0], [0.942809041, 0, 0.3333333], [-0.47140452, 0.81649658, 0.33333333],
                               [-0.47140452, -0.81649658, 0.33333333], [0, 0, -1]]),
                             ("prolate_rotors", "CH3Cl",
                              [[0, 0, 0.68], [0.942809042, 0, 1.01333334], [-0.47140452, 0.81649658, 1.01333334],
                               [-0.47140452, -0.81649658, 1.01333334], [0, 0, -0.319999995]]),
                             ("asymmetric_rotors", "H2O",
                              [[-0.02618914, -0.04536092, -0.03703704], [0.9166199, -0.04536092, 0.2962963],
                               [-0.49759366, 0.77113566, 0.2962963]])
                         ])
def test_flip_up_down(rotor_type, atomic_protein, atom_coord, request):
    atomic_prot_obj = request.getfixturevalue(rotor_type)[atomic_protein]
    atomic_prot_obj.flip_up_down()
    assert np.allclose([atom.coord for atom in atomic_prot_obj.pdb.get_atoms()], atom_coord, atol=precision)


@pytest.mark.parametrize("rotor_type, atomic_protein",
                         [("spherical_rotors", "CH4")]
                         )
def test_flatten_error(rotor_type, atomic_protein, request):
    atomic_prot_obj = request.getfixturevalue(rotor_type)[atomic_protein]
    with pytest.raises(ValueError) as excinfo:
        atomic_prot_obj.flatten()
    print("\nError:", str(excinfo.value))


@pytest.mark.parametrize("rotor_type, atomic_protein, com, I, atom_coord",
                         [
                             ("linear_rotors", "H2", "H2_com", "H2_I", [[0, 0, 0.8660254], [0, 0, -0.8660254]]),
                             ("linear_rotors", "HD", "HD_com", "HD_I", [[0, 0, 1.15470053], [0, 0, -0.577350281]]),
                             ("prolate_rotors", "CH3Cl", "CH3Cl_com", "CH3Cl_I",
                              [[0, 0, -0.68], [np.sqrt(8 / 9), 0, -1 / 3 - 0.68],
                               [-np.sqrt(2 / 9), np.sqrt(2 / 3), -1 / 3 - 0.68],
                               [-np.sqrt(2 / 9), - np.sqrt(2 / 3), -1 / 3 - 0.68],
                               [0, 0, 1 - 0.68]
                               ]),
                             ("asymmetric_rotors", "H2O", "H2O_com", "H2O_I",
                              [[-0.0582641549, -0.0268424014, 0],
                               [0.46611324, 0.21473922, 0.81649658], [0.46611324, 0.21473922, -0.81649658]])
                         ])
def test_flatten(rotor_type, atomic_protein, com, I, atom_coord, request):
    atomic_prot_obj = request.getfixturevalue(rotor_type)[atomic_protein]
    atomic_prot_obj.flatten()
    assert np.allclose([atom.coord for atom in atomic_prot_obj.pdb.get_atoms()], atom_coord, atol=precision)


@pytest.mark.parametrize("rotor_type, atomic_protein",
                         [
                             ("asymmetric_rotors", "H2O")
                         ])
def test_skew_sym_matrix(rotor_type, atomic_protein, request):
    fixture = request.getfixturevalue(rotor_type)
    test = [[
        [0, 0, 0], np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])],
        [[1, 1, 1], np.array([[0, -1, 1], [1, 0, -1], [-1, 1, 0]])]
    ]
    for row in test:
        assert np.allclose(row[1], fixture[atomic_protein]._skew_sym_matrix(row[0]), atol=precision)


@pytest.mark.parametrize("rotor_type, atomic_protein, mass",
                         [("asymmetric_rotors", "H2O", "H2O_mass"),
                          ])
def test_mass(rotor_type, atomic_protein, mass, request):
    fixture = request.getfixturevalue(rotor_type)
    assert np.isclose(fixture[atomic_protein].mass, fixture[mass], atol=precision)


@pytest.mark.parametrize("rotor_type, atomic_protein, name",
                         [("asymmetric_rotors", "H2O", "H2O_name"),
                          ])
def test_name(rotor_type, atomic_protein, name, request):
    fixture = request.getfixturevalue(rotor_type)
    assert fixture[atomic_protein].name == fixture[name]


@pytest.mark.parametrize("rotor_type, atomic_protein, rot, tran, result_com, result_coord",
                         [
                             ("linear_rotors", "H2", np.eye(3), [1, 2, 3], [1, 2, 3],
                              [[0.5, 1.5, 2.5], [1.5, 2.5, 3.5]]),  # transl v = [1,2,3]
                             ("linear_rotors", "H2", np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]),
                              np.zeros(3), [0, 0, 0], [[-0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),  # Rx 90deg
                             ("linear_rotors", "H2", np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]),
                              np.zeros(3), [0, 0, 0], [[0.5, -0.5, -0.5], [-0.5, 0.5, 0.5]]),  # Ry 90deg
                             ("linear_rotors", "H2", np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]),
                              np.zeros(3), [0, 0, 0], [[0.5, 0.5, -0.5], [-0.5, -0.5, 0.5]]),  # Ry 180,
                             ("linear_rotors", "H2", np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]), [1, 2, 3], [1, 2, 3],
                              [[0.5, 1.5, 3.5], [1.5, 2.5, 2.5]]),  # Rx 90deg + transl v = [1,2,3]
                             ("asymmetric_rotors", "H2O_shifted", np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]),
                              [1, 2, 3], [1, 2, 3], [[0.9738108, 2.0370370, 3.0453609],
                                                     [1.916619, 1.703703, 3.0453609],
                                                     [0.5024063, 1.703703, 2.2288643]])  # Rx 90deg + transl v = [1,2,3]
                         ])
def test_transform_pdb(rotor_type, atomic_protein, rot, tran, result_com, result_coord, request):
    fixture = request.getfixturevalue(rotor_type)
    ap = fixture[atomic_protein]
    ap.transform_pdb(rot=rot, tran=tran)
    assert np.allclose([atom.coord for atom in ap.pdb.get_atoms()], result_coord, atol=precision)


def test_prepare_pigments():
    def _build_atoms():
        atom0 = FakeAtom((-1, -1, 0), mass=1, element='N', id=0)  # NA
        atom1 = FakeAtom((1, -1, 0), mass=1, element='N', id=1)  # ND
        atom2 = FakeAtom((-1, 1, 0), mass=1, element='N', id=2)  # NB
        atom3 = FakeAtom((1, 1, 0), mass=1, element='N', id=3)  # NC
        atom4 = FakeAtom((0, 0, 0), mass=1, element='MG', id=4)  # all masses eq to 1 here

        return {"NA": atom0, "ND": atom1, "NB": atom2, "NC": atom3, "MG": atom4}

    def _build_FakeResidue():
        CLA = FakeResidue(dict_atoms=_build_atoms(), name="CLA")
        CLB = FakeResidue(dict_atoms=_build_atoms(), name="CLB")
        CLC = FakeResidue(dict_atoms=_build_atoms(), name="CLC")
        CLD = FakeResidue(dict_atoms=_build_atoms(), name="CLD")

        return {"CLA": CLA, "CLB": CLB, "CLC": CLC, "CLD": CLD}

    dict_residues = _build_FakeResidue()
    ppc_pdb = FakePDB([atom for res in dict_residues.values() for atom in res.get_atoms()], id=1,
                      child_list=[res for res in dict_residues.values()])

    ap = ProteinAtomic(ppc_pdb, name='fake_pigment_protein')
    ap_2 = ProteinAtomic(ppc_pdb, name='fake_pigment_protein2')
    dict_fake_residues = _build_FakeResidue()
    for resname in dict_fake_residues.keys():
        ap.prepare_pigments(resname, ChlorophyllAtomic)
        ap_2.prepare_pigments(resname, ChlorophyllAtomic, 'data_for_testing')

    tempDict = {}
    for key, value in dict_fake_residues.items():
        tempDict[f'{value.get_parent().id}_{key}_{value.id[1]}'] = ChlorophyllAtomic(
            value)

    for key, value in tempDict.items():
        assert ap.dict_pigments[key].residue.name == value.residue.name
        assert ap.dict_pigments[key].residue.id[1] == value.residue.id[1]
        assert ap_2.dict_pigments[key].residue.name == value.residue.name
        assert ap_2.dict_pigments[key].residue.id[1] == value.residue.id[1]
        assert ap_2.dict_pigments[key].dict_data == {'tresp_atoms': ['NA', 'NB', 'MG'], 'tresp_pc': [0.25, 0.25, 0.5],
                                                     'vacuum_mag': 0.8, 'dipole_mag': 1.5}


def test_construct_subcomplex():
    # prepare pdb objects with one model
    H2_xyz_pdb = FakePDB([FakeAtom([-0.5, -0.5, -0.5]), FakeAtom([0.5, 0.5, 0.5])], id=1, list_chains=['X', 'Y', 'Z'])
    H2_xy_pdb = FakePDB([FakeAtom([-0.5, -0.5, -0.5]), FakeAtom([0.5, 0.5, 0.5])], id=1, list_chains=['X', 'Y'])
    H2_y_pdb = FakePDB([FakeAtom([-0.5, -0.5, -0.5]), FakeAtom([0.5, 0.5, 0.5])], id=1, list_chains=['Y'])
    H2_atomic = ProteinAtomic(pdb=H2_xyz_pdb, name="H2_full")
    subset_atomic_one = H2_atomic.construct_subcomplex(['X', 'Y'], 'XY_chains')
    subset_atomic_two = H2_atomic.construct_subcomplex(['Y'], 'X_chain')
    # tests
    assert f'{subset_atomic_one.pdb.child_list[0].child_list}' == f'{H2_xy_pdb.child_list[0].child_list}'
    assert f'{subset_atomic_two.pdb.child_list[0].child_list}' == f'{H2_y_pdb.child_list[0].child_list}'
    assert f'{H2_atomic.pdb.child_list[0].child_list}' == f'{H2_xyz_pdb.child_list[0].child_list}'
    with pytest.raises(ValueError) as excinfo:
        H2_atomic.construct_subcomplex(['A'], 'Bad_chain_call')
    print("\nError:", str(excinfo.value))

    # prepare pdb objects with multiple models
    H2_pdb = FakePDB([FakeAtom([-0.5, -0.5, -0.5]), FakeAtom([0.5, 0.5, 0.5])], id=1, list_chains=['X', 'Y', 'Z'],
                     child_list=[FakeChild(['X', 'Y'], id=0), FakeChild(['X', 'Y'], id=1)])
    H2_pdb_x1 = FakePDB([FakeAtom([-0.5, -0.5, -0.5]), FakeAtom([0.5, 0.5, 0.5])], id=1, list_chains=['X', 'Y', 'Z'],
                        child_list=[FakeChild(['X', 'Y'], id=0), FakeChild(['X'], id=1)])
    H2_twomodel_atomic = ProteinAtomic(H2_pdb, 'two_models')
    H2_twomodel_sub = H2_twomodel_atomic.construct_subcomplex(['X'], 'X_chain', 1)
    # tests
    assert f'{H2_twomodel_sub.pdb.child_list[0].child_list}' == f'{H2_pdb_x1.child_list[0].child_list}'
    with pytest.raises(ValueError) as excinfo:
        H2_twomodel_atomic.construct_subcomplex(['X'], 'X_chain')
    print("\nError:", str(excinfo.value))


@pytest.mark.parametrize("rotor_type, atomic_protein, name",
                         [("asymmetric_rotors", "H2O", "H2O_name")])
def test_name(rotor_type, atomic_protein, name, request):
    fixture = request.getfixturevalue(rotor_type)
    assert fixture[atomic_protein].name == fixture[name]


def test_load_hamiltonian():
    def _build_atoms():
        atom0 = FakeAtom((-1, -1, 0), mass=1, element='N', id=0)  # NA
        atom1 = FakeAtom((1, -1, 0), mass=1, element='N', id=1)  # ND
        atom2 = FakeAtom((-1, 1, 0), mass=1, element='N', id=2)  # NB
        atom3 = FakeAtom((1, 1, 0), mass=1, element='N', id=3)  # NC
        atom4 = FakeAtom((0, 0, 0), mass=1, element='MG', id=4)  # all masses eq to 1 here

        return {"NA": atom0, "ND": atom1, "NB": atom2, "NC": atom3, "MG": atom4}

    def _build_FakeResidue():
        CLA = FakeResidue(dict_atoms=_build_atoms(), name="CLA")
        CLB = FakeResidue(dict_atoms=_build_atoms(), name="CLB")
        CLC = FakeResidue(dict_atoms=_build_atoms(), name="CLC")

        return {"CLA": CLA, "CLB": CLB, "CLC": CLC}

    dict_residues = _build_FakeResidue()
    CH4_pdb = FakePDB([atom for res in dict_residues.values() for atom in res.get_atoms()], id=1,
                      child_list=[res for res in dict_residues.values()])

    # object_0
    CH4_Atomic_0 = ProteinAtomic(CH4_pdb, name='fake_pigment_protein')
    for resname in dict_residues.keys():
        CH4_Atomic_0.prepare_pigments(resname, ChlorophyllAtomic)
    # object_1
    CH4_Atomic_1 = ProteinAtomic(CH4_pdb, name='fake_pigment_protein')
    for resname in dict_residues.keys():
        CH4_Atomic_1.prepare_pigments(resname, ChlorophyllAtomic)

    # load fake Hamiltonian:
    path = os.path.abspath(__file__)[:-(len('testing/test_atomic_protein.py'))]
    H_fake_same_pdb_label = pd.read_csv(path + 'util/hamiltonian_data/test_H_same_PDB_label.csv')
    H_fake_diff_pdb_label = pd.read_csv(path + 'util/hamiltonian_data/test_H_different_PDB_label.csv')
    # A) In the case of dict_pigments_map_H_PDB, dict_map_chains=None
    # ================================================================
    CH4_Atomic_0.load_hamiltonian('testing_dict_h_same_PDB_label', dict_pigments_map_H_PDB=None,
                                  dict_map_chains=None)

    # check if the H label are loaded correctly:
    assert set(CH4_Atomic_0.dict_data['hamiltonian'].columns.values) == set(H_fake_same_pdb_label.columns.values)
    # check if the H values are loaded correctly:
    assert np.array_equal(CH4_Atomic_0.dict_data['hamiltonian'].values, H_fake_same_pdb_label.values)

    with pytest.raises(ValueError) as excinfo:
        CH4_Atomic_0.load_hamiltonian('testing_dict_h_diff_PDB_label')
    print("\nError in the case of dict_pigments_map_H_PDB is None: \n", str(excinfo.value))

    # # B) In the case of dict_pigments_map_H_PDB, dict_map_chains is not None:
    # # # =================================================================
    CH4_Atomic_1.load_hamiltonian('testing_dict_h_diff_PDB_label',
                                  dict_pigments_map_H_PDB={"XXX": "ID1_CLA_()", "YYY": "ID1_CLB_()",
                                                           "ZZZ": "ID1_CLC_()"},
                                  dict_map_chains={"ID1": "X"})
    # check if the H label are loaded correctly:
    assert set(CH4_Atomic_1.dict_data['hamiltonian'].columns.values) == set(['X_CLA_()', 'X_CLB_()', 'X_CLC_()'])
    # check if the H values are loaded correctly:
    assert np.array_equal(CH4_Atomic_1.dict_data['hamiltonian'].values, H_fake_diff_pdb_label.values)

    # test mapping cases:
    # ===================
    # dict_pigments_map_H_PDB is not complete:

    with pytest.raises(ValueError) as excinfo:
        CH4_Atomic_1.load_hamiltonian('testing_dict_h_diff_PDB_label',
                                      dict_pigments_map_H_PDB={"XXX": "ID1_CLA_()", "ZZZ": "ID1_CLC_()"},
                                      dict_map_chains={"ID1": "y", 'ID2': 'd'})
    print("\nError in the case of dict_pigments_map_H_PDB is Not None: \n", str(excinfo.value))

    # if dict_pigments_map_H_PDB was constructed wrongly:
    with pytest.raises(ValueError) as excinfo:
        CH4_Atomic_1.load_hamiltonian('testing_dict_h_diff_PDB_label',
                                      dict_pigments_map_H_PDB={"XXX": "ID2_CLA_()", "YYY": "ID1_CLB_()",
                                                               "ZZZ": "ID1_CLC_()"},
                                      dict_map_chains={"ID1": "y"})
    print("\nError in the case of dict_pigments_map_H_PDB is Not None: \n", str(excinfo.value))
