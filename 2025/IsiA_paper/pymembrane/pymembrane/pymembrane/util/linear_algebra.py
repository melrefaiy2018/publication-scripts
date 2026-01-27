import numpy as np
from scipy.spatial.transform import Rotation as R

from pymembrane.util.physical_constants import precision


def construct_rotation_matrix(angle, axis):

    if axis not in ['x', 'y', 'z']:
        raise ValueError('Axis must be x, y, or z.')

    r = R.from_euler(axis, angle, degrees=False)
    return r.as_matrix()


def skew_sym_matrix(V1_vec):
    """
    This helper function calculates the skew-symmetric cross product that we need
    for rotating the 'special axis' of the moment of inertia into the z-axis.

    Parameters
    ----------
    1. V1_vec : np.array
                the input array

    Returns
    -------
    1. V2_skew_sym_cross : np.array
                           the skew sym cross array
    """
    return np.array([[0, -V1_vec[2], V1_vec[1]],
                     [V1_vec[2], 0, -V1_vec[0]],
                     [-V1_vec[1], V1_vec[0], 0]])


def construct_rotation_for_flatten(I2_inertia, Z1_unit=None):
    if Z1_unit is None:
        Z1_unit = np.array([0, 0, 1])
    else:
        Z1_unit = np.array(Z1_unit)
        if not np.isclose(np.linalg.norm(Z1_unit), 1):
            raise ValueError('Z1_unit must be a unit vector.')


    # Make inertia zero if really close to zero so that eigen value is valid
    I2_inertia[np.isclose(I2_inertia, 0, atol=precision)] = 0

    w, v = np.linalg.eig(I2_inertia)
    # For spherical rotors, return value error because you can not flatten a spherical rotor
    if np.all(np.isclose(w, w[0])):
        raise ValueError("spherical rotors cannot do flatten")
    w_dist = [np.abs(w[0] - w[1]), np.abs(w[0] - w[2]), np.abs(w[2] - w[1])]
    if np.argmin(w_dist) == 0:
        index_z = 2
    elif np.argmin(w_dist) == 1:
        index_z = 1
    else:
        index_z = 0

    # Now rotate the vector v[index_z] into the z-axis
    # The following procedure is non-obvious and described here:
    # http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    V1_z = v[:, index_z]
    #  if v[index_z] already on z-axis skip the rotate and transform it back
    if np.all(V1_z == Z1_unit):
        return np.eye(3)
    V1_cross = np.cross(V1_z, Z1_unit)
    V2_x = skew_sym_matrix(V1_cross)
    R2_rot = np.identity(3) + V2_x + np.dot(V2_x, V2_x) * (
                1 - np.dot(V1_z, Z1_unit)) / np.linalg.norm(
        V1_cross) ** 2

    return R2_rot
