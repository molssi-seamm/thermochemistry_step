# -*- coding: utf-8 -*-

"""Thermochemistry keeps each species' Hessian (Paul 2026-11-26):
hessians/<SMILES>.hessian.npz in the step's directory and the property
hessian file#Thermochemistry#<model>."""

import types

import numpy as np
from ase import Atoms
from ase.vibrations import VibrationsData

from seamm_exec.hessian import frequencies, load_hessian
from seamm_util import Q_
from thermochemistry_step.cbh_approach import CBHApproach


class Properties(dict):
    def exists(self, name):
        return name in self

    def add(self, name, _type, units=None, description=None):
        self[name] = {"value": None}

    def put(self, name, value):
        self[name] = {"value": value}

    def get(self, name):
        return self[name]


def test_the_species_hessian_is_kept(tmp_path):
    root = tmp_path / "job"
    step = root / "3"
    step.mkdir(parents=True)
    atoms = Atoms("H2", positions=[[0, 0, 0], [0.74, 0, 0]])
    k = 0.5  # E_h/bohr^2 along the bond
    H = np.zeros((6, 6))
    H[0, 0] = H[3, 3] = k
    H[0, 3] = H[3, 0] = -k
    H_ev = Q_(H, "hartree/bohr**2").m_as("eV/Å**2")
    me = types.SimpleNamespace(
        _vibrations=VibrationsData(atoms, H_ev.reshape(2, 3, 2, 3)),
        _species_directory=lambda: step,
        flowchart=types.SimpleNamespace(root_directory=root),
        model="ORCA:CMP@R2SCAN-3C",
        _fd_info=None,
        _hessian_check={"errors": [0.001, 0.002, 0.001]},
    )
    properties = Properties()
    properties.put("frequency fingerprint#Thermochemistry#ORCA:CMP@R2SCAN-3C", "fp")
    configuration = types.SimpleNamespace(
        properties=properties, charge=0, spin_multiplicity=1
    )
    relative = CBHApproach._save_species_hessian(
        me, configuration, "[H][H]", "analytic"
    )
    assert relative == "3/hessians/[H][H].hessian.npz"
    assert properties["hessian file#Thermochemistry#ORCA:CMP@R2SCAN-3C"] == {
        "value": relative
    }
    data = load_hessian(root / relative)
    assert np.allclose(data["hessian"], H)
    assert data["symbols"] == ["H", "H"]
    info = data["info"]
    assert info["route"] == "analytic" and info["settings"] == "fp"
    assert info["species"] == "[H][H]" and info["check"]["errors"][1] == 0.002
    # the frequencies from the file are the species' own
    nu = frequencies(data["hessian"], data["masses"])
    assert abs(nu.max() - max(me._vibrations.get_frequencies().real)) < 0.5


def test_keeping_never_fails_the_species(tmp_path):
    me = types.SimpleNamespace(
        _vibrations=object(),  # not a VibrationsData
        _species_directory=lambda: tmp_path,
        flowchart=None,
        model="m",
    )
    configuration = types.SimpleNamespace(
        properties=Properties(), charge=0, spin_multiplicity=1
    )
    assert CBHApproach._save_species_hessian(me, configuration, "C", "x") is None
