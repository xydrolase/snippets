#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
geometry.py
Utility module implementing 3d geometrical operations related to
the internal and external coordinates of polymer structure.

Author: Xin Yin <xinyin at iastate dot edu>
"""

__version__ = 0.1


import math

import numpy as np
import numpy.matlib

def rodrigues_rot_matrix(axis, phi):
    """Use Rodrigues' formula to compute a rotation matrix corresponding
    to a rotation by angle phi about a fixed axis specified by vector
    w = <x, y, z>."""

    def anti_symmetric_matrix(x, y, z):
        return np.matrix([
            [0, -z, y],
            [z, 0, -x],
            [-y, x, 0]
            ])

    w = np.array(axis) / np.linalg.norm(axis)
    Omega = anti_symmetric_matrix(*list(w))

    return np.matlib.eye(3) + math.sin(phi) * Omega + \
            2 * math.sin(phi/2) ** 2 * (Omega * Omega)

def rotate_about(v, rot_axis, phi):
    if type(v) in (list, tuple):
        v = np.matrix(v).transpose()
    if type(v) is np.ndarray and len(v.shape) == 1:
        v = np.matrix(v).transpose()
    else:
        raise TypeError("Invalid dimension of v.")

    rot_v = rodrigues_rot_matrix(rot_axis, phi) * v

    return np.squeeze(np.asarray(rot_v))

def normal_p3(p1, p2, p3):
    return normal_v2(p1-p2, p3-p2)

def normal_v2(v1, v2):
    n = np.cross(v1, v2)
    return n / np.linalg.norm(n)
