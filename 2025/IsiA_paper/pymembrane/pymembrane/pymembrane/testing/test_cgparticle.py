import numpy as np
import pytest
from testing.test_atomic_pigment import FakeResidue
from testing.test_atomic_protein import FakeAtom, FakePDB

from pymembrane.structure import cg_particle, atomic_protein, atomic_pigment
from pymembrane.structure.atomic_pigment import ChlorophyllAtomic
from pymembrane.structure.atomic_protein import ProteinAtomic
from pymembrane.structure.cg_particle import CGProtein
from pymembrane.util.physical_constants import precision

"""
CGProtein
"""


def rotation_matrices():
    """basic rotations around x, y, z axes"""
    phi = np.pi / 2
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)

    """Rotation matrices  in 3D"""
    Rx = np.array([[1, 0, 0],
                   [0, cos_phi, -sin_phi],
                   [0, sin_phi, cos_phi]])

    Ry = np.array([[cos_phi, 0, sin_phi],
                   [0, 1, 0],
                   [-sin_phi, 0, cos_phi]])

    Rz = np.array([[cos_phi, -sin_phi, 0],
                   [sin_phi, cos_phi, 0],
                   [0, 0, 1]])
    Rxyz = Rz @ Ry @ Rx  # all the three rotation in sequence

    return {"Rx": Rx, "Ry": Ry, "Rz": Rz, "Rxyz": Rxyz}


def translation_matrices():
    v_transl = np.array((1, 1, 1))
    w_transl = np.array((-1, -1, -1))
    return {"v_transl": v_transl, "w_transl": w_transl}


def positions():
    p000 = (0, 0, 0)
    p111 = (1, 1, 1)
    return {"p000": p000, "p111": p111}


def _build_atoms():
    atom0 = FakeAtom((-1, -1, 0), mass=1, element='N', id=0)  # NA
    atom1 = FakeAtom((1, -1, 0), mass=1, element='N', id=1)  # ND
    atom2 = FakeAtom((-1, 1, 0), mass=1, element='N', id=2)  # NB
    atom3 = FakeAtom((1, 1, 0), mass=1, element='N', id=3)  # NC
    atom4 = FakeAtom((0, 0, 0), mass=1, element='MG', id=4)  # all masses eq to 1 here

    return {"NA": atom0, "ND": atom1, "NB": atom2, "NC": atom3, "MG": atom4}


def _build_protein_atomic():
    protein_pdb = FakePDB([_build_atoms()['NA'], _build_atoms()['ND'], _build_atoms()['NB'],
                           _build_atoms()['NC'], _build_atoms()['MG']])  # this contains 1 residues
    protein_atomic = atomic_protein.ProteinAtomic(pdb=protein_pdb, name='name_atomic_protein')
    return protein_atomic


def _build_protein_atomic_with_pigments():
    """
    Build protein atomic that includes pigments
    """
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
        return {"CLA": CLA, "CLB": CLB}
    dict_residues = _build_FakeResidue()
    CH4_pdb = FakePDB([atom for res in dict_residues.values() for atom in res.get_atoms()], id=1,
                      child_list=[res for res in dict_residues.values()])
    protein_Atomic = ProteinAtomic(CH4_pdb, name='fake_pigment_protein')
    for resname in dict_residues.keys():
        protein_Atomic.prepare_pigments(resname, ChlorophyllAtomic)

    return protein_Atomic



def _build_protein_cg():
    return {
        "p000_i": cg_particle.CGProtein(_build_protein_atomic(), "p000_i", location=positions()["p000"],
                                        R2_orient=np.eye(3)),
        "p111_i": cg_particle.CGProtein(_build_protein_atomic(), "p111_i", location=positions()["p111"],
                                        R2_orient=np.eye(3)),
        "p000_rx": cg_particle.CGProtein(_build_protein_atomic(), "p000_rx", location=positions()["p000"],
                                         R2_orient=rotation_matrices()["Rx"]),
        "p000_ry": cg_particle.CGProtein(_build_protein_atomic(), "p000_ry", location=positions()["p000"],
                                         R2_orient=rotation_matrices()["Ry"]),
        "p000_rz": cg_particle.CGProtein(_build_protein_atomic(), "p000_rz", location=positions()["p000"],
                                         R2_orient=rotation_matrices()["Rz"]),
        "p000_rxyz": cg_particle.CGProtein(_build_protein_atomic(), "p000_rxyz", location=positions()["p000"],
                                           R2_orient=rotation_matrices()["Rxyz"]),
        "p111_rx": cg_particle.CGProtein(_build_protein_atomic(), "p111_rx", location=positions()["p111"],
                                         R2_orient=rotation_matrices()["Rx"]),
        "p111_ry": cg_particle.CGProtein(_build_protein_atomic(), "p111_ry", location=positions()["p111"],
                                         R2_orient=rotation_matrices()["Ry"]),
        "p111_rz": cg_particle.CGProtein(_build_protein_atomic(), "p111_rz", location=positions()["p111"],
                                         R2_orient=rotation_matrices()["Rz"]),
        "p111_rxyz": cg_particle.CGProtein(_build_protein_atomic(), "p111_rxyz", location=positions()["p111"],
                                           R2_orient=rotation_matrices()["Rxyz"]),
    }


def _protein_cg_translated(protein_cg, v):
    protein_cg.translate(v)
    return protein_cg


def _protein_cg_rotated_by_matrix(protein_cg, R):
    protein_cg.rotate_by_matrix(R)
    return protein_cg


def _protein_cg_rotated_by_axis(protein_cg, angle, axis):
    protein_cg.rotate_by_axis(angle, axis)
    return protein_cg




@pytest.mark.parametrize(
    "computed, expected",
    [(_build_protein_cg()["p000_i"].location, positions()["p000"]),
     (_build_protein_cg()["p111_i"].location, positions()["p111"]),
     (_build_protein_cg()["p000_i"].R2_orient, np.eye(3)),
     (_build_protein_cg()["p000_rx"].R2_orient, rotation_matrices()["Rx"]),
     (_build_protein_cg()["p000_ry"].R2_orient, rotation_matrices()["Ry"]),
     (_build_protein_cg()["p000_rz"].R2_orient, rotation_matrices()["Rz"]),
     ])
def test_init_CGProtein(computed, expected):
    assert np.allclose(computed, expected, atol=precision)



@pytest.mark.parametrize(
    "computed, expected",
    [(_build_protein_cg()["p000_i"].name, "p000_i")])
def test_name_CGProtein(computed, expected):
    assert computed == expected


@pytest.mark.parametrize(
    'computed, expected',
    [(_protein_cg_translated(_build_protein_cg()["p000_i"], translation_matrices()["v_transl"]).location, (1, 1, 1)),
     (_protein_cg_translated(_build_protein_cg()["p000_i"], translation_matrices()["w_transl"]).location, (-1, -1, -1)),
     (_protein_cg_translated(_build_protein_cg()["p111_i"], translation_matrices()["v_transl"]).location, (2, 2, 2)),
     (_protein_cg_translated(_build_protein_cg()["p000_rx"], translation_matrices()["v_transl"]).location, (1, 1, 1)),
     (_protein_cg_translated(_build_protein_cg()["p000_ry"], translation_matrices()["v_transl"]).location, (1, 1, 1)),
     (_protein_cg_translated(_build_protein_cg()["p000_rz"], translation_matrices()["v_transl"]).location, (1, 1, 1))])
def test_translate_GCProtein(computed, expected):
    assert np.allclose(computed, expected, atol=precision)


@pytest.mark.parametrize(
    'computed, expected',
    [(_protein_cg_translated(_build_protein_cg()["p000_i"], translation_matrices()["v_transl"]).location, (1, 1, 1)),
     (_protein_cg_translated(_build_protein_cg()["p000_i"], translation_matrices()["w_transl"]).location, (-1, -1, -1)),
     (_protein_cg_translated(_build_protein_cg()["p111_i"], translation_matrices()["v_transl"]).location, (2, 2, 2)),
     (_protein_cg_translated(_build_protein_cg()["p000_rx"], translation_matrices()["v_transl"]).location, (1, 1, 1)),
     (_protein_cg_translated(_build_protein_cg()["p000_ry"], translation_matrices()["v_transl"]).location, (1, 1, 1)),
     (_protein_cg_translated(_build_protein_cg()["p000_rz"], translation_matrices()["v_transl"]).location, (1, 1, 1))])
def test_translate_GCProtein(computed, expected):
    assert np.allclose(computed, expected, atol=precision)


@pytest.mark.parametrize(
    'computed, expected',
    [(_protein_cg_rotated_by_matrix(_build_protein_cg()["p000_i"], np.eye(3)).R2_orient, np.eye(3)),
     (_protein_cg_rotated_by_matrix(_build_protein_cg()["p111_i"], np.eye(3)).R2_orient, np.eye(3)),
     (_protein_cg_rotated_by_matrix(_build_protein_cg()["p000_rx"], np.eye(3)).R2_orient, rotation_matrices()["Rx"]),
     (_protein_cg_rotated_by_matrix(_build_protein_cg()["p000_ry"], np.eye(3)).R2_orient, rotation_matrices()["Ry"]),
     (_protein_cg_rotated_by_matrix(_build_protein_cg()["p000_rz"], np.eye(3)).R2_orient, rotation_matrices()["Rz"]),
     (_protein_cg_rotated_by_matrix(_build_protein_cg()["p000_rxyz"], np.eye(3)).R2_orient,
      rotation_matrices()["Rxyz"])])
def test_rotation_by_matrix_GCProtein(computed, expected):
    assert np.allclose(computed, expected, atol=precision)


@pytest.mark.parametrize(
    'computed, expected',
    [(_protein_cg_rotated_by_axis(_build_protein_cg()["p000_i"], 0, 'x').R2_orient, np.eye(3)),
     (_protein_cg_rotated_by_axis(_build_protein_cg()["p111_i"], 0, 'x').R2_orient, np.eye(3)),
     (_protein_cg_rotated_by_axis(_build_protein_cg()["p000_i"], 0, 'y').R2_orient, np.eye(3)),
     (_protein_cg_rotated_by_axis(_build_protein_cg()["p000_i"], 0, 'z').R2_orient, np.eye(3)),
     (_protein_cg_rotated_by_axis(_build_protein_cg()["p000_i"], np.pi / 2, 'x').R2_orient, rotation_matrices()["Rx"]),
     (_protein_cg_rotated_by_axis(_build_protein_cg()["p111_i"], np.pi / 2, 'y').R2_orient, rotation_matrices()["Ry"]),
     (_protein_cg_rotated_by_axis(_build_protein_cg()["p111_i"], np.pi / 2, 'z').R2_orient, rotation_matrices()["Rz"])])
def test_rotation_by_axis_GCProtein(expected, computed):
    with pytest.raises(ValueError):
        _protein_cg_rotated_by_axis(_build_protein_cg()["p000_i"], 0, 'randomstring')
    assert np.allclose(expected, computed, atol=precision)


def test_atomic_protein_in_init():
    p0_atomic = _build_protein_atomic_with_pigments()
    protein_cg = CGProtein(p0_atomic, name='protein_CG1')
    assert protein_cg.atomic is p0_atomic

def test_name_in_init():
    p0_atomic = _build_protein_atomic_with_pigments()
    protein_cg = CGProtein(p0_atomic, name='protein_CG1')
    assert protein_cg.name == 'protein_CG1'


def test_dict_pigment_in_init():
    p0_atomic = _build_protein_atomic_with_pigments()
    protein_cg = CGProtein(p0_atomic, name='protein_CG1')
    assert protein_cg.dict_pigments.keys() == p0_atomic.dict_pigments.keys()
    assert [pig.atomic for pig in protein_cg.dict_pigments.values()] == list(
        p0_atomic.dict_pigments.values())

########################################################################################################################

"""
CGPigment
"""


def _build_pigment_atomic():
    return atomic_pigment.ChlorophyllAtomic(FakeResidue(_build_atoms()))


def _build_pigment_cg():
    return {
        "p000_i": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p000_i"]),
        "p111_i": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p111_i"]),
        "p000_rx": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p000_rx"]),
        "p000_ry": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p000_ry"]),
        "p000_rz": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p000_rz"]),
        "p000_rxyz": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p000_rxyz"]),
        "p111_rx": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p111_rx"]),
        "p111_ry": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p111_ry"]),
        "p111_rz": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p111_rz"]),
        "p111_rxyz": cg_particle.CGPigment(_build_pigment_atomic(), _build_protein_cg()["p111_rxyz"]),
    }

def test_init_CGPigment():
    p0_atomic = _build_protein_atomic_with_pigments()
    protein_cg = CGProtein(p0_atomic, name='protein_CG1')
    cg_pig = list(protein_cg.dict_pigments.values())[0]
    atomic_pig = list(p0_atomic.dict_pigments.values())[0]
    assert cg_pig.cg_protein is protein_cg
    assert cg_pig.atomic is atomic_pig


@pytest.mark.parametrize("computed, expected",
                         [(_build_pigment_cg()["p000_i"].get_dipole_dir(), [-0.70710678, 0.70710678, 0]),
                          (_build_pigment_cg()["p000_rx"].get_dipole_dir(),
                           [-0.707106781, 4.32978028e-17, 0.707106781]),
                          (_build_pigment_cg()["p000_ry"].get_dipole_dir(),
                           [-4.32978028e-17, 0.707106781, 0.707106781]),
                          (_build_pigment_cg()["p000_rz"].get_dipole_dir(), [-0.70710678, -0.70710678, 0])
                          ])
def test_get_dipole_dir(computed, expected):
    assert np.allclose(computed, expected, atol=precision)


@pytest.mark.parametrize(
    "computed, expected",
    [(_build_pigment_cg()["p000_i"].mass, _build_pigment_atomic().mass),
     (_build_pigment_cg()["p111_i"].mass, _build_pigment_atomic().mass),
     (_build_pigment_cg()["p000_rx"].mass, _build_pigment_atomic().mass),
     (_build_pigment_cg()["p000_ry"].mass, _build_pigment_atomic().mass)])
def test_mass_CGPigment(computed, expected):
    assert np.allclose(computed, expected, atol=precision)




@pytest.mark.parametrize(
    "computed, expected",
    [(_build_pigment_cg()["p000_i"].xyz, [0, 0, 0]),
     (_build_pigment_cg()["p111_i"].xyz, [1, 1, 1]),
     (_build_pigment_cg()["p000_rx"].xyz, [0, 0, 0]),
     (_build_pigment_cg()["p000_ry"].xyz, [0, 0, 0]),
     (_build_pigment_cg()["p000_rz"].xyz, [0, 0, 0]),
     (_build_pigment_cg()["p000_rxyz"].xyz, [0, 0, 0]),
     (_build_pigment_cg()["p111_rx"].xyz, [1, 1, 1]),
     (_build_pigment_cg()["p111_ry"].xyz, [1, 1, 1])
     ])
def test_xyz_CGPigment(computed, expected):
    assert np.allclose(computed, expected, atol=precision)


@pytest.mark.parametrize(
    "computed, expected",
    [(_build_pigment_cg()["p000_i"].get_atom_xyz('NA'), [-1, -1, 0]),
     (_build_pigment_cg()["p111_i"].get_atom_xyz('NA'), [0, 0, 1]),
     (_build_pigment_cg()["p000_rx"].get_atom_xyz('NA'), [-1, 0, -1]),
     (_build_pigment_cg()["p000_ry"].get_atom_xyz('NA'), [0, -1, 1]),
     (_build_pigment_cg()["p000_rz"].get_atom_xyz('NA'), [1, -1, 0]),
     (_build_pigment_cg()["p111_rx"].get_atom_xyz('NA'), [0, 1, 0]),
     (_build_pigment_cg()["p111_ry"].get_atom_xyz('NA'), [1, 0, 2])
     ])
def test_get_atom_xyz(computed, expected):

    assert np.allclose(computed, expected, atol=precision)
