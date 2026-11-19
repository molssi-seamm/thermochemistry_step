# -*- coding: utf-8 -*-

"""The lowest few force-field conformers of a species are optimized at the
geometry level and the lowest there is used (design 2026-11-19): UFF puts
gauche diarsine 0.2 kJ/mol below trans, but trans is the ground state (CCSD(T)
1.4 kJ/mol lower, r2SCAN-3c 2.0), and the step used to build gauche. The
geometry level is a stub here: trans is 2.0 kJ/mol lower."""

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from seamm_util import Q_  # noqa: E402

DIARSINE = "[AsH2][AsH2]"


def _P(**changes):
    return {
        "initial conformer": "lowest of the embeddings (MMFF)",
        "target conformer": "as given",
        "conformer embeddings": 20,
        "conformers optimized": 3,
        "optimize fragments": "yes",
        "thermal corrections": "harmonic",
        "fragment convergence": Q_(0.005, "eV/Å"),
        "fragment max steps": 100,
        **changes,
    }


def _is_trans(configuration):
    """Whether diarsine is trans: two H-As-As-H dihedrals of 180 degrees."""
    symbols = list(configuration.atoms.symbols)
    xyz = np.array(configuration.atoms.get_coordinates(fractionals=False))
    a, b = [i for i, s in enumerate(symbols) if s == "As"]
    hs = [i for i, s in enumerate(symbols) if s == "H"]
    ha = [h for h in hs if np.linalg.norm(xyz[h] - xyz[a]) < 1.8]
    hb = [h for h in hs if h not in ha]

    def dihedral(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 = b1 / np.linalg.norm(b1)
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        return abs(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))

    angles = [dihedral(xyz[x], xyz[a], xyz[b], xyz[y]) for x in ha for y in hb]
    return sum(1 for d in angles if d > 170.0) == 2


def _energy(configuration):
    return -2.0 if _is_trans(configuration) else 0.0


def _diarsine(node):
    db = node._test_db
    system = db.create_system(name=f"CBH:{DIARSINE}", make_current=True)
    configuration = system.create_configuration(name="initial", make_current=True)
    configuration.from_smiles(DIARSINE)
    configuration.charge = 0
    configuration.spin_multiplicity = 1
    return configuration


def test_embed_lowest_gives_the_candidates(node):
    configuration = _diarsine(node)
    candidates = node._embed_lowest(configuration, DIARSINE, 20, keep=3)
    assert len(candidates) == 2 and candidates[0][0] == 0.0
    assert configuration.n_atoms == 6 and configuration.charge == 0
    # the configuration is at the first, the force field's lowest: gauche
    first = np.array(configuration.atoms.get_coordinates(fractionals=False))
    assert np.allclose(first, candidates[0][1])
    assert not _is_trans(configuration)
    # one candidate only: nothing to rank
    assert node._embed_lowest(configuration, DIARSINE, 20) is None
    assert node._embed_lowest(configuration, "C", 20, keep=3) is None


def test_the_lowest_at_the_geometry_level_is_used(node, monkeypatch):
    """Through the subflowchart: _optimize leaves its energy."""
    configuration = _diarsine(node)
    candidates = node._embed_lowest(configuration, DIARSINE, 20, keep=3)
    seen = []

    def optimize(configuration, _P, can="the species", fmax=None):
        node._last_energy = _energy(configuration)
        seen.append(node._last_energy)
        return True

    monkeypatch.setattr(node, "_optimize", optimize)
    chosen = node._lowest_at_geometry_level(_P(), DIARSINE, configuration, candidates)
    assert chosen == 1 and seen == [0.0, -2.0]
    assert _is_trans(configuration)


def test_the_model_chemistry_ranks_them_too(node, monkeypatch):
    configuration = _diarsine(node)
    candidates = node._embed_lowest(configuration, DIARSINE, 20, keep=3)
    node._engine_mode = "Model Chemistry"

    def optimize_and_energy(configuration, optimize=True, fmax=0.005, max_steps=200):
        node._last_optimization = {"converged": True}
        return _energy(configuration)

    monkeypatch.setattr(node, "optimize_and_energy", optimize_and_energy)
    assert node._lowest_at_geometry_level(_P(), DIARSINE, configuration, candidates)
    assert _is_trans(configuration)


def test_without_an_engine_the_force_field_lowest_is_kept(node, monkeypatch):
    """With the Model Chemistry run as separate tasks nothing is optimized, and
    energies of force-field structures would rank nothing."""
    configuration = _diarsine(node)
    candidates = node._embed_lowest(configuration, DIARSINE, 20, keep=3)
    node._engine_mode = "Model Chemistry"
    calls = []

    def energy_only(configuration, optimize=True, fmax=0.005, max_steps=200):
        calls.append(1)
        return _energy(configuration)

    monkeypatch.setattr(node, "optimize_and_energy", energy_only)
    assert (
        node._lowest_at_geometry_level(_P(), DIARSINE, configuration, candidates) == 0
    )
    assert calls == [1] and not _is_trans(configuration)


def test_how_many_are_optimized(node):
    assert node._conformers_to_optimize(_P()) == 3
    assert node._conformers_to_optimize(_P(**{"conformers optimized": 1})) == 1
    assert node._conformers_to_optimize(_P(**{"optimize fragments": "no"})) == 1
    # an older flowchart, without the parameter
    P = _P()
    del P["conformers optimized"]
    assert node._conformers_to_optimize(P) == 3
    # the subflowchart's energy without harmonic corrections is not optimized
    no_harmonic = _P(**{"thermal corrections": "none (the energy is the enthalpy)"})
    assert node._conformers_to_optimize(no_harmonic) == 1
    node._engine_mode = "Model Chemistry"
    assert node._conformers_to_optimize(no_harmonic) == 3


def test_a_fragment_is_computed_in_its_lowest_conformer(node, monkeypatch):
    """_compute_species builds the fragment from its SMILES, ranks the
    candidates, and computes the lowest."""
    node._engine_mode = "Model Chemistry"

    def optimize_and_energy(configuration, optimize=True, fmax=0.005, max_steps=200):
        node._last_optimization = {"converged": True}
        return _energy(configuration)

    computed = []

    def compute(_P, can, configuration):
        computed.append(_is_trans(configuration))
        return {"E": _energy(configuration)}

    monkeypatch.setattr(node, "optimize_and_energy", optimize_and_energy)
    monkeypatch.setattr(node, "_compute_species_mc", compute)
    monkeypatch.setattr(node, "_check_species", lambda *a: None)
    monkeypatch.setattr(node, "_write_species_sdf", lambda *a: None)
    result = node._compute_species(_P(), DIARSINE, DIARSINE)
    assert computed == [True] and result["E"] == -2.0

    # and with one conformer optimized, the force field's (gauche)
    computed.clear()
    monkeypatch.setattr(node, "_find_fragment", lambda db, can: None)
    node._compute_species(_P(**{"conformers optimized": 1}), DIARSINE, DIARSINE)
    assert computed == [False]


def test_the_table_shows_rotamers_of_hydrogens(node):
    """The candidates' dihedrals include hydrogens, so trans and gauche diarsine,
    and anti and gauche ethanol, can be told apart (thermo, 5837)."""
    configuration = _diarsine(node)
    candidates = node._embed_lowest(configuration, DIARSINE, 20, keep=3)
    symbols = list(configuration.atoms.symbols)
    index = {_id: k for k, _id in enumerate(configuration.atoms.ids)}
    bonds = [(index[r["i"]], index[r["j"]]) for r in configuration.bonds.bonds()]
    texts = [node._rotamer_dihedrals(symbols, bonds, xyz) for _, xyz in candidates]
    assert any(t.startswith("As1-As2 ") and t.endswith("/180/180") for t in texts)
    assert any(not t.endswith("/180/180") for t in texts)
    # Ethanol: only C-C-O-H; the methyl rotor and the ring-free C-H bonds are left out
    symbols = ["C", "C", "O", "H", "H", "H", "H", "H", "H"]
    bonds = [(0, 1), (1, 2), (0, 3), (0, 4), (0, 5), (1, 6), (1, 7), (2, 8)]
    xyz = [
        [0.0, 0.0, 0.0],
        [1.52, 0.0, 0.0],
        [2.0, 1.35, 0.0],
        [-0.4, 1.0, 0.0],
        [-0.4, -0.5, 0.87],
        [-0.4, -0.5, -0.87],
        [1.9, -0.5, 0.87],
        [1.9, -0.5, -0.87],
        [2.97, 1.35, 0.0],
    ]
    assert node._rotamer_dihedrals(symbols, bonds, xyz) == "C2-O3 180"
    # and a ring bond is no rotor
    assert node._rotamer_dihedrals(["C"] * 3, [(0, 1), (1, 2), (0, 2)], xyz[:3]) == ""
