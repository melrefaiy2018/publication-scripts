import pytest
from testing.test_atomic_protein import FakeAtom, FakeResidue

from pymembrane.structure.atomic_pigment import *
from pymembrane.util.physical_constants import precision


@pytest.fixture(scope="module")
def build_atoms():
    # building a 5 atoms fake pigment
    # ===============================
    atom0 = FakeAtom((0, 0, 0), mass=14, element='N', id=0)  # NA
    atom1 = FakeAtom((0, 2, 0), mass=14, element='N', id=1)  # ND
    atom2 = FakeAtom((2, 0, 0), mass=14, element='N', id=2)  # NB
    atom3 = FakeAtom((2, 2, 0), mass=14, element='N', id=3)  # NC
    atom4 = FakeAtom((1, 1, 0), mass=24, element='MG', id=4)
    dict_atoms = {'NA': atom0, 'NB': atom2, 'NC': atom3, 'ND': atom1, 'MG': atom4}
    residue_5aa = FakeResidue(dict_atoms, name="test_name")
    chl = ChlorophyllAtomic(residue_5aa)
    chl_mass = 80

    chl_xyz = np.array((1, 1, 0))  # unitary
    chl_qx = np.array([-2, -2, 0]) / np.linalg.norm(np.array([-2, -2, 0]))
    chl_qy = np.array([2, -2, 0]) / np.linalg.norm(np.array([2, -2, 0]))

    # defining bond based on different scaling factor (s.f.):
    # ================================================
    chl_list_of_bonds_0 = ((0, 4), (1, 4), (2, 4), (3, 4))  # default s.f. = 1.05
    chl_list_of_bonds_1 = ((0, 1), (2, 4), (1, 2), (0, 4), (3, 4), (0, 3), (1, 4),
                           (2, 3), (0, 2), (1, 3))  # s.f. 2.7
    chl_list_of_bonds_2 = (())  # s.f. 0.6

    yield {"chl": chl, "chl_mass": chl_mass, "chl_xyz": chl_xyz, "chl_qx": chl_qx,
           "chl_qy": chl_qy, "atom0": atom0, "atom1": atom1, "atom2": atom2,
           "atom3": atom3, "atom4": atom4, "chl_list_of_bonds_0": chl_list_of_bonds_0,
           "chl_list_of_bonds_1": chl_list_of_bonds_1,
           "chl_list_of_bonds_2": chl_list_of_bonds_2, "residue": residue_5aa}

@pytest.fixture(scope="module")
def build_pheophytin():
    # building a 5 atoms fake pigment
    # ===============================
    atom0 = FakeAtom((0, 0, 0), mass=14, element='N', id=0)  # NA
    atom1 = FakeAtom((0, 2, 0), mass=14, element='N', id=1)  # ND
    atom2 = FakeAtom((2, 0, 0), mass=14, element='N', id=2)  # NB
    atom3 = FakeAtom((2, 2, 0), mass=14, element='N', id=3)  # NC
    dict_atoms = {'NA': atom0, 'NB': atom2, 'NC': atom3, 'ND': atom1,}
    residue_4aa = FakeResidue(dict_atoms, name="test_pheophytin")
    PHO = PheophytinAtomic(residue_4aa)
    pho_mass = 80

    pho_xyz = np.array((1, 1, 0))  # unitary

    # defining bond based on different scaling factor (s.f.):
    # ================================================
    pho_list_of_bonds_0 = ((0, 4), (1, 4), (2, 4), (3, 4))  # default s.f. = 1.05


    yield {"pho": PHO, "pho_mass": pho_mass, "pho_xyz": pho_xyz,
           "atom0": atom0, "atom1": atom1, "atom2": atom2,
           "atom3": atom3, "pho_list_of_bonds_0": pho_list_of_bonds_0,
           "residue": residue_4aa}


@pytest.mark.parametrize("function, atomic_pigment, dipole, a1, a2, res, bonds",
                         [
                             ("build_atoms", "chl", "qx", "atom0", "atom3",
                              "residue", 'chl_list_of_bonds_0'),
                             ("build_atoms", "chl", "qy", "atom2", "atom1",
                              'residue', 'chl_list_of_bonds_0')
                         ])
def test_initialization(function, atomic_pigment, dipole, a1, a2, res, bonds, request):
    fixture = request.getfixturevalue(function)
    if dipole == "qx":
        assert (fixture[atomic_pigment].qx_atoms == (fixture[a1], fixture[a2]))
    if dipole == "qy":
        assert (fixture[atomic_pigment].qy_atoms == (fixture[a1], fixture[a2]))
    assert sorted(fixture[atomic_pigment].bonds) == sorted(fixture[bonds])
    assert fixture[atomic_pigment].dict_data == {}
    assert fixture[atomic_pigment].residue == fixture[res]

@pytest.mark.parametrize("function, atomic_pigment, mass",
                         [
                             ("build_atoms", "chl", "chl_mass")
                         ])
def test_mass(function, atomic_pigment, mass, request):
    fixture = request.getfixturevalue(function)
    assert np.isclose(fixture[atomic_pigment].mass, fixture[mass], atol=precision)


@pytest.mark.parametrize("function, atomic_pigment, xyz",
                         [
                             ("build_atoms", "chl", "chl_xyz"),
                             ("build_pheophytin", "pho", "pho_xyz")
                         ])
def test_xyz(function, atomic_pigment, xyz, request):
    fixture = request.getfixturevalue(function)
    assert np.allclose(fixture[atomic_pigment].xyz, fixture[xyz], atol=precision)


@pytest.mark.parametrize("function, atomic_pigment, param, expected",
                         [
                             ("build_atoms", "chl", None, "chl_qy"),
                             ("build_atoms", "chl", "qy", "chl_qy"),
                             ("build_atoms", "chl", "qx", "chl_qx"),
                             ("build_atoms", "chl", "randomstring", None)
                         ])
def test_get_dipole_dir(function, atomic_pigment, param, expected, request):
    fixture = request.getfixturevalue(function)
    if param is None:
        assert np.allclose(fixture[atomic_pigment].get_dipole_dir(), fixture[expected], atol=precision)
    elif expected is None and (param != "qx" and param != "qy"):
        with pytest.raises(ValueError):
            fixture[atomic_pigment].get_dipole_dir(param)
    else:
        assert np.allclose(fixture[atomic_pigment].get_dipole_dir(param), fixture[expected], atol=precision)


@pytest.mark.parametrize("function, atomic_pigment, scale_factor, list_bonds",
                         [
                             ("build_atoms", "chl", None, "chl_list_of_bonds_0"),
                             ("build_atoms", "chl", "2.7", "chl_list_of_bonds_1"),
                             ("build_atoms", "chl", "0.6", "chl_list_of_bonds_2")
                         ])
def test_chlorophyll_define_bonds(function, atomic_pigment, scale_factor, list_bonds, request):
    fixture = request.getfixturevalue(function)
    if scale_factor is None:
        assert np.allclose(sorted(fixture[atomic_pigment].bonds), sorted(fixture[list_bonds]), atol=precision)
    else:
        assert np.allclose(sorted(fixture[atomic_pigment].define_bonds(scale_factor=float(scale_factor))),
                            sorted(fixture[list_bonds]), atol=precision)


@pytest.mark.parametrize("function, atomic_pigment, expected",
                         [
                             ("build_atoms", "chl", "X_test_name_()")
                         ])
def test_name(function, atomic_pigment, expected, request):
    fixture = request.getfixturevalue(function)
    assert fixture[atomic_pigment].name == expected


@pytest.mark.parametrize("function, atomic_pigment,atom_str, expected",
                         [
                             ("build_atoms", "chl", "NA", [0, 0, 0]),
                             ("build_atoms", "chl", "NB", [2, 0, 0]),
                             ("build_atoms", "chl", "NC", [2, 2, 0]),
                             ("build_atoms", "chl", "ND", [0, 2, 0]),
                             ("build_atoms", "chl", "MG", [1, 1, 0]),
                         ])
def test_get_atom_xyz(function, atomic_pigment, atom_str, expected, request):
    fixture = request.getfixturevalue(function)
    assert np.allclose(fixture[atomic_pigment].get_atom_xyz(atom_str), expected, atol=precision)


def test_get_coupling_data():
    atom0 = FakeAtom((0, 0, 0), mass=14, element='N', id=0)  # NA
    atom1 = FakeAtom((0, 2, 0), mass=14, element='N', id=1)  # ND
    atom2 = FakeAtom((2, 0, 0), mass=14, element='N', id=2)  # NB
    atom3 = FakeAtom((2, 2, 0), mass=14, element='N', id=3)  # NC
    atom4 = FakeAtom((1, 1, 0), mass=24, element='MG', id=4)
    dict_atoms = {'NA': atom0, 'NB': atom2, 'NC': atom3, 'ND': atom1, 'MG': atom4}
    residue_5aa = FakeResidue(dict_atoms, name="test_name")
    residue_5aa.dict_data = {'tresp_atoms': ['NA', 'NB', 'MG'], 'tresp_pc': [0.25, 0.25, 0.5],
                             'vacuum_mag': 0.8, 'dipole_mag': 1.5}
    chl = ChlorophyllAtomic(residue_5aa)
    chl.get_coupling_data('data_for_testing')
    assert residue_5aa.dict_data == chl.dict_data
    with pytest.raises(KeyError) as excinfo:
        chl.get_coupling_data('bad_dictionary_name')
    print("\nError:", str(excinfo.value))

