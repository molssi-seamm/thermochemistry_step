# -*- coding: utf-8 -*-

"""Reusing a target stored in the fragment store, only when the stored structure
is this target's: the same atoms in the same order, the same bonds, and the
heavy atoms within 0.3 Å RMSD with no reflection (Paul, 2026-11-13: "careful
that the structures do indeed match")."""

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from thermochemistry_step.cbh_approach import kabsch_rmsd, match_structure  # noqa

from tests.test_cbh_energy_level import (  # noqa: E402,F401
    ENERGIES,
    FakeEvaluator,
    LOW,
    OFFSET,
    _computed,
    split,
)

#: Butane, anti (heavy-atom dihedral 180), and Cartesian helpers
BUTANE = ["C", "C", "C", "C"] + ["H"] * 10


def _rotate(xyz, axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array(
        [[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]]
    )
    r = np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k
    return np.asarray(xyz) @ r.T


def _butane(dihedral=180.0):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    mol = Chem.AddHs(Chem.MolFromSmiles("CCCC"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, 0, 1, 2, 3, dihedral)
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    return symbols, np.array(conf.GetPositions())


def test_kabsch_ignores_rotation_and_translation():
    a = np.random.default_rng(1).normal(size=(6, 3))
    b = _rotate(a, [1, 2, 3], 1.1) + [5.0, -2.0, 1.0]
    assert kabsch_rmsd(a, b) == pytest.approx(0.0, abs=1e-10)


def test_the_same_structure_matches():
    symbols, xyz = _butane()
    moved = _rotate(xyz, [0, 1, 1], 0.7) + 3.0
    noise = np.random.default_rng(2).normal(scale=0.02, size=xyz.shape)
    ok, rmsd, why = match_structure(symbols, xyz, symbols, moved + noise)
    assert ok and rmsd < 0.06 and why is None


def test_another_conformer_does_not_match():
    symbols, anti = _butane(180.0)
    _, gauche = _butane(60.0)
    ok, rmsd, why = match_structure(symbols, anti, symbols, gauche)
    assert not ok and rmsd > 0.3 and "another conformer" in why


def test_another_order_or_other_bonds_do_not_match():
    symbols, xyz = _butane()
    order = list(range(len(symbols)))
    order[0], order[3] = order[3], order[0]
    swapped = [symbols[i] for i in order]
    assert swapped == symbols  # the same symbols, other atoms at each index
    ok, _, why = match_structure(symbols, xyz, swapped, xyz[order])
    assert not ok and why == "other bonds"
    ok, _, why = match_structure(symbols, xyz, symbols[::-1], xyz[::-1])
    assert not ok and why == "other atoms or atom order"
    broken = xyz.copy()
    broken[3] += [3.0, 0.0, 0.0]  # pull a terminal carbon off
    ok, _, why = match_structure(symbols, xyz, symbols, broken)
    assert not ok and why == "other bonds"


def test_the_mirror_image_does_not_match():
    """CHFClBr and its enantiomer: the same bonds, but no proper rotation brings
    one onto the other."""
    symbols = ["C", "H", "F", "Cl", "Br"]
    xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.63, 0.63, 0.63],
            [-0.78, -0.78, 0.78],
            [-1.02, 1.02, -1.02],
            [1.11, -1.11, -1.11],
        ]
    )
    mirror = xyz * [1, 1, -1]
    ok, rmsd, _ = match_structure(symbols, xyz, symbols, mirror)
    assert not ok and rmsd > 0.3
    ok, _, _ = match_structure(symbols, xyz, symbols, xyz)
    assert ok


def _store_target(path, split, shift=0.0):  # noqa: F811
    """Store the target CCF's geometry as given (moved by `shift` Å along x for
    one hydrogen-free change: a translation, which must still match)."""
    from seamm_thermochemistry import EnergyRecord, FragmentStore, GeometryRecord

    _, configuration = split.get_system_configuration()
    symbols = list(configuration.atoms.symbols)
    xyz = np.asarray(configuration.atoms.get_coordinates(fractionals=False)) + shift
    E = ENERGIES["CCF"] + OFFSET["CCF"]
    with FragmentStore(path, create=True) as store:
        g = GeometryRecord.make(
            "CCF",
            symbols,
            xyz.tolist(),
            E=E,
            ZPE=10.0,
            H=E + 15.0,
            model_chemistry=LOW,
            frequency_fingerprint="fp:low",
            optimization_fingerprint="fp:low",
            version="6.1.1",
        )
        store.put_geometry(g)
        store.put_energy(
            EnergyRecord.make(
                g["geometry_hash"], ENERGIES["CCF"], energy_fingerprint="fp:high"
            )
        )
    return symbols, xyz


def _predicts(node, monkeypatch, versions):
    """The program predicts 'fp:low' at version 6.1.x, records the versions asked."""

    def expected(structure, mc, kind, version=None):
        versions.append(version)
        return ["fp:low"] if str(version or "").startswith("6.1") else None

    monkeypatch.setattr(node, "_expected_fingerprints", expected)


def test_a_stored_target_is_reused(split, tmp_path, monkeypatch):  # noqa: F811
    path = tmp_path / "fragments.db"
    _store_target(path, split, shift=1.5)
    versions = []
    _predicts(split, monkeypatch, versions)
    P = _computed(split, **{"fragment store": "use and add what is computed"})
    P["fragment store file"] = str(path)
    data = split.run_cbh(P)
    assert "CCF" not in split._computed
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    # The same result as computing it
    P["reuse stored target"] = "never"
    split._computed.clear()
    again = split.run_cbh(P)
    assert "CCF" in split._computed
    assert again["DfH 298K"] == pytest.approx(data["DfH 298K"], abs=1e-6)


def test_another_stored_structure_is_not_reused(split, tmp_path):  # noqa: F811
    """A stored target whose structure is not this one's (here: its atoms in
    another order) is computed."""
    from seamm_thermochemistry import FragmentStore, GeometryRecord

    path = tmp_path / "fragments.db"
    _, configuration = split.get_system_configuration()
    symbols = list(configuration.atoms.symbols)
    xyz = np.asarray(configuration.atoms.get_coordinates(fractionals=False))
    order = list(range(len(symbols)))[::-1]
    E = ENERGIES["CCF"] + OFFSET["CCF"]
    with FragmentStore(path, create=True) as store:
        store.put_geometry(
            GeometryRecord.make(
                "CCF",
                [symbols[i] for i in order],
                xyz[order].tolist(),
                E=E,
                ZPE=10.0,
                H=E + 15.0,
                model_chemistry=LOW,
            )
        )
    P = _computed(split, **{"fragment store": "use and add what is computed"})
    P["fragment store file"] = str(path)
    split.run_cbh(P)
    assert "CCF" in split._computed


def test_unpredictable_settings_are_not_reused(split, tmp_path):  # noqa: F811
    """Without the program's prediction of its fingerprints the settings of the
    stored target cannot be checked, so it is computed."""
    path = tmp_path / "fragments.db"
    _store_target(path, split)
    P = _computed(split, **{"fragment store": "use and add what is computed"})
    P["fragment store file"] = str(path)
    split.run_cbh(P)
    assert "CCF" in split._computed


def test_another_version_fails_the_target(split, tmp_path, monkeypatch):  # noqa: F811
    """A target taken from the store at one version of the program, in a run
    whose calculations report another, fails instead of mixing versions."""
    split._assumed_version = "6.1.1"
    split._fingerprints = {"energy": {"fingerprint": "fp:x", "version": "6.2.0"}}
    _, configuration = split.get_system_configuration()
    with pytest.raises(RuntimeError, match="computed with version 6.1.1"):
        split._store_fingerprints(configuration, None, "CCF")
    split._fingerprints = {"energy": {"fingerprint": "fp:x", "version": "6.1.2"}}
    split._store_fingerprints(configuration, None, "CCF")  # the same 6.1
