# -*- coding: utf-8 -*-

"""The rotational symmetry number of a molecule, perceived from its geometry.

The symmetry number σ is the number of proper rotations (including the
identity) in the molecule's point group: the number of indistinguishable
orientations reached by rotation alone. It enters the rotational partition
function as 1/σ, i.e. the entropy as -R ln σ. This module finds σ by testing
candidate rotation axes directly against the structure, which needs no
point-group assignment and no external library.

Candidate axes are the principal axes of inertia, the vectors from the centre
of mass to each atom, the midpoints and cross products of pairs of
like atoms, and the normals of planes through triples of like atoms. For every
candidate axis and every order n from 2 up to the largest number of atoms of
one element, the rotation by 2π/n is applied and accepted if it maps every
atom onto an atom of the same element within the tolerance.
"""

from collections import Counter
from itertools import combinations

import numpy as np

__all__ = ["rotational_symmetry_number"]


def _rotation_matrix(axis, angle):
    """Rodrigues' formula for the rotation by `angle` about the unit vector `axis`."""
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    C = 1.0 - c
    return np.array(
        [
            [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
        ]
    )


def _maps_onto_itself(xyz, groups, R, tol):
    """True if the rotation R maps every atom onto a like atom within tol."""
    rotated = xyz @ R.T
    for indices in groups.values():
        a = rotated[indices]
        b = xyz[indices]
        # Each rotated atom must have a partner; with equal counts and a
        # rotation (a bijection on the molecule) one partner each suffices.
        d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        if np.any(d.min(axis=1) > tol):
            return False
    return True


def _unique_axes(candidates, tol=1.0e-3):
    """Normalize the candidate vectors and drop (anti)parallel duplicates."""
    axes = []
    for v in candidates:
        norm = np.linalg.norm(v)
        if norm < 1.0e-6:
            continue
        u = v / norm
        if any(abs(np.dot(u, w)) > 1.0 - tol for w in axes):
            continue
        axes.append(u)
    return axes


def rotational_symmetry_number(symbols, xyz, masses=None, tol=0.1, max_order=12):
    """The rotational symmetry number of a molecule.

    Parameters
    ----------
    symbols : [str]
        The element symbol of each atom.
    xyz : array-like, shape (n_atoms, 3)
        The Cartesian coordinates, Å.
    masses : array-like, optional
        The atomic masses, for the centre of mass. Equal masses are used if
        omitted (the centre of mass then coincides with the centroid; for
        the symmetry this only matters through the candidate axes).
    tol : float
        The distance, in Å, within which a rotated atom must land on a like
        atom.
    max_order : int
        The largest rotation order tested.

    Returns
    -------
    sigma : int
        The symmetry number.
    axes : {int: int}
        The number of distinct proper rotation axes of each order found,
        e.g. ``{6: 1, 2: 7}`` for benzene (the C6 axis also carries C3 and
        C2, which are not counted separately), ``{3: 4, 2: 3}`` for methane.
    """
    xyz = np.asarray(xyz, dtype=float)
    n_atoms = len(symbols)
    if n_atoms == 1:
        return 1, {}

    if masses is None:
        masses = np.ones(n_atoms)
    masses = np.asarray(masses, dtype=float)
    xyz = xyz - (masses[:, None] * xyz).sum(axis=0) / masses.sum()

    groups = {}
    for i, symbol in enumerate(symbols):
        groups.setdefault(symbol, []).append(i)
    groups = {k: np.array(v) for k, v in groups.items()}

    # Linear molecules: sigma is 2 if inversion through the centre of mass
    # is a symmetry operation, else 1. Collinearity: all cross products with
    # the longest atom vector vanish.
    lengths = np.linalg.norm(xyz, axis=1)
    longest = xyz[np.argmax(lengths)]
    if np.all(np.linalg.norm(np.cross(xyz, longest), axis=1) < tol * lengths.max()):
        if _maps_onto_itself(xyz, groups, -np.eye(3), tol):
            return 2, {"linear": 2}
        return 1, {"linear": 1}

    # Candidate axes
    candidates = []
    inertia = np.zeros((3, 3))
    for m, r in zip(masses, xyz):
        inertia += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    _, vectors = np.linalg.eigh(inertia)
    candidates.extend(vectors.T)
    candidates.extend(xyz)
    for indices in groups.values():
        if len(indices) < 2:
            continue
        for i, j in combinations(indices, 2):
            candidates.append(xyz[i] + xyz[j])
            candidates.append(np.cross(xyz[i], xyz[j]))
        # Normals of planes through three like atoms (a C_n axis with n >= 3
        # is normal to the plane of each orbit of n atoms). For a populous
        # element only triples of near neighbours are used, since the atoms
        # of an orbit are adjacent in its polygon.
        if len(indices) <= 20:
            for i, j, k in combinations(indices, 3):
                candidates.append(np.cross(xyz[j] - xyz[i], xyz[k] - xyz[i]))
        else:
            sub = xyz[indices]
            d = np.linalg.norm(sub[:, None, :] - sub[None, :, :], axis=2)
            for a in range(len(indices)):
                nearest = np.argsort(d[a])[1:7]
                for b, c in combinations(nearest, 2):
                    candidates.append(np.cross(sub[b] - sub[a], sub[c] - sub[a]))
    axes = _unique_axes(candidates)

    highest = min(max_order, max(Counter(symbols).values()))

    operations = set()
    found = {}
    for axis in axes:
        best = 1
        for n in range(2, highest + 1):
            R = _rotation_matrix(axis, 2.0 * np.pi / n)
            if _maps_onto_itself(xyz, groups, R, tol):
                best = max(best, n)
                for k in range(1, n):
                    Rk = _rotation_matrix(axis, 2.0 * np.pi * k / n)
                    operations.add(tuple(np.round(Rk, 3).ravel()))
        if best > 1:
            found[best] = found.get(best, 0) + 1

    sigma = 1 + len(operations)
    return sigma, dict(sorted(found.items(), reverse=True))
