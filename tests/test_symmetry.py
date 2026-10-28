#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the perceived rotational symmetry number."""

import numpy as np
import pytest
from ase.build import molecule

from thermochemistry_step.symmetry import rotational_symmetry_number


@pytest.mark.parametrize(
    "name, sigma",
    [
        ("H2O", 2),
        ("NH3", 3),
        ("CH4", 12),
        ("C6H6", 12),
        ("C2H6", 6),
        ("C2H4", 4),
        ("CO2", 2),
        ("HCl", 1),
        ("CH3Cl", 3),
        ("C2H2", 2),
        ("H2", 2),
        ("CO", 1),
        ("H2O2", 2),
        ("CH3OH", 1),
        ("BF3", 6),
        ("PH3", 3),
        ("SiH4", 12),
        ("N2O", 1),
        ("H2CO", 2),
        ("CH3CHO", 1),
        ("C60", 60),
    ],
)
def test_g2_molecules(name, sigma):
    atoms = molecule(name)
    result, _ = rotational_symmetry_number(
        atoms.get_chemical_symbols(), atoms.positions, atoms.get_masses()
    )
    assert result == sigma


def test_octahedral_SF6():
    xyz = np.array(
        [
            [0, 0, 0],
            [1.6, 0, 0],
            [-1.6, 0, 0],
            [0, 1.6, 0],
            [0, -1.6, 0],
            [0, 0, 1.6],
            [0, 0, -1.6],
        ]
    )
    sigma, axes = rotational_symmetry_number(["S"] + ["F"] * 6, xyz)
    assert sigma == 24
    assert axes == {4: 3, 3: 4, 2: 6}


def test_single_atom():
    assert rotational_symmetry_number(["Ar"], [[0.0, 0.0, 0.0]]) == (1, {})


def test_slightly_distorted_geometry_within_tolerance():
    atoms = molecule("C6H6")
    rng = np.random.default_rng(1)
    positions = atoms.positions + rng.normal(0.0, 0.005, atoms.positions.shape)
    sigma, _ = rotational_symmetry_number(
        atoms.get_chemical_symbols(), positions, atoms.get_masses()
    )
    assert sigma == 12


def test_distorted_geometry_beyond_tolerance_loses_symmetry():
    atoms = molecule("C6H6")
    rng = np.random.default_rng(1)
    positions = atoms.positions + rng.normal(0.0, 0.2, atoms.positions.shape)
    sigma, _ = rotational_symmetry_number(
        atoms.get_chemical_symbols(), positions, atoms.get_masses(), tol=0.02
    )
    assert sigma == 1
