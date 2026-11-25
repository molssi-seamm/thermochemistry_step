# -*- coding: utf-8 -*-

"""A broken Hessian is not a saddle point (Thermochemistry 6043, 2026-11-25):
4-hydroxy-TEMPO's analytic Hessian had a frequency of -3.1e7 cm^-1 and was chased
as a saddle for four hours. A Hessian with any frequency above MAX_FREQUENCY is
computed again once, at the same structure in a new engine session; if it is
still broken, a target fails and a fragment is flagged (its rungs fail, it is not
stored), with no displacement."""

import contextlib

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from seamm_util import Q_  # noqa: E402

CH4 = [
    ("C", (0.0, 0.0, 0.0)),
    ("H", (0.629, 0.629, 0.629)),
    ("H", (-0.629, -0.629, 0.629)),
    ("H", (-0.629, 0.629, -0.629)),
    ("H", (0.629, -0.629, -0.629)),
]
#: A Hessian (E_h/bohr^2) whose frequencies are all below MAX_FREQUENCY
GOOD = np.diag(np.linspace(0.05, 0.4, 15))
#: The same with one element absurdly negative (a frequency of ~1e5i cm^-1)
BAD = GOOD.copy()
BAD[3, 3] = -5.0e3
BOHR = 0.529177210903
MC = {
    "level": "ORCA:DFT@B3LYP/def2-SVP",
    "step": "ORCA",
    "method": "B3LYP",
    "options": {"mdi_capable": True, "analytic_hessian": True},
}


class Engine:
    def __init__(self, hessians):
        self.x = np.zeros(15)
        self.hessians = list(hessians)
        self.n_hessians = 0

    def set_coordinates(self, x, units="Å"):
        self.x = np.asarray(x, dtype=float).reshape(-1) / BOHR

    def energy(self, units="kJ/mol"):
        e = 0.5 * self.x @ GOOD @ self.x
        return e * 2625.4996394799 if units == "kJ/mol" else e

    def supports(self, command):
        return command == "<HESSIAN"

    def hessian(self):
        self.n_hessians += 1
        return self.hessians.pop(0).copy()

    def fingerprints(self):
        return None


@pytest.fixture()
def methane(tmp_path):
    from molsystem import SystemDB

    def make(name):
        db = SystemDB(filename=str(tmp_path / f"methane-{name.replace(':', '_')}.db"))
        configuration = db.create_system(name=name).create_configuration()
        configuration.atoms.append(
            x=[p[0] for _, p in CH4],
            y=[p[1] for _, p in CH4],
            z=[p[2] for _, p in CH4],
            symbol=[s for s, _ in CH4],
        )
        return configuration

    return make


def _run(node, configuration, monkeypatch, hessians):
    import seamm_exec.hessian as H

    engine = Engine(hessians)
    sessions = []

    @contextlib.contextmanager
    def session(self, configuration, announce=True):
        xyz = np.array(configuration.atoms.get_coordinates(fractionals=False))
        s = H.EngineSession(self, engine, configuration, None, [6, 1, 1, 1, 1], xyz)

        def optimize(fmax=0.005, max_steps=200):
            self.last_optimization = {"converged": True, "evaluations": 1}
            return True

        s.optimize = optimize
        sessions.append(s)
        yield s

    def thermochemistry(_P, *args, **kwargs):
        return {
            "total_energy": Q_(node._E0, "eV").m_as("E_h"),
            "ZPE": 100.0,
            "H thermal": [[110.0]],
        }

    displaced = []
    monkeypatch.setattr(H.HessianDriver, "session", session)
    monkeypatch.setattr(node, "model_chemistry", lambda: MC)
    monkeypatch.setattr(node, "harmonic_thermochemistry", thermochemistry)
    monkeypatch.setattr(
        node,
        "_largest_imaginary_mode",
        lambda *a, **k: displaced.append(1) or (None, None),
    )
    node._hessian_route = "warm engine over MDI"
    node.model = MC["level"]
    P = {
        "thermal corrections": "harmonic",
        "optimize fragments": "yes",
        "fragment convergence": Q_(0.005, "eV/Å"),
        "fragment max steps": 200,
        "step size": Q_(0.01, "Å"),
    }
    result = node._compute_species_mc(P, "C", configuration)
    return result, engine, sessions, displaced


def test_a_sound_hessian_is_used_as_is(node, methane, monkeypatch):
    result, engine, sessions, _ = _run(node, methane("C"), monkeypatch, [GOOD])
    assert engine.n_hessians == 1 and len(sessions) == 1
    assert "hessian refused" not in result and result["ZPE"] == 100.0


def test_a_broken_hessian_is_computed_again(node, methane, monkeypatch):
    result, engine, sessions, _ = _run(node, methane("C"), monkeypatch, [BAD, GOOD])
    # once more, at the same structure, in a new engine session
    assert engine.n_hessians == 2 and len(sessions) == 2
    assert "hessian refused" not in result and result["ZPE"] == 100.0


def test_a_target_whose_hessian_stays_broken_fails(node, methane, monkeypatch):
    with pytest.raises(RuntimeError, match="Hessian of C is broken"):
        _run(node, methane("C"), monkeypatch, [BAD, BAD])


def test_a_fragment_whose_hessian_stays_broken_is_flagged(node, methane, monkeypatch):
    result, engine, sessions, displaced = _run(
        node, methane("CBH:C"), monkeypatch, [BAD, BAD]
    )
    assert engine.n_hessians == 2
    assert "frequency of" in result["hessian refused"]
    assert "ZPE" not in result and "H298" not in result
    # no saddle handling on a broken Hessian
    assert displaced == []
    # the rungs that need it fail, and it is not stored
    from thermochemistry_step.cbh_approach import CBHApproach

    class R:
        coefficients = {"C": 1, "CC": -1}

    failed = CBHApproach._failed_rungs("CC", [1], {1: R()}, {"C": result, "CC": {}})
    assert failed == {1: ["C"]}
    assert CBHApproach._storable("C", result) is False


def test_the_refusal_names_the_frequency(node, methane):
    from ase import Atoms
    from ase.vibrations import VibrationsData

    atoms = Atoms("CH4", positions=[p for _, p in CH4])
    for K, refused in ((GOOD, False), (BAD, True)):
        h = Q_(K, "hartree/bohr**2").m_as("eV/Å**2")
        node._vibrations = VibrationsData(atoms, h.reshape(5, 3, 5, 3))
        reason = node._hessian_refusal()
        assert (reason is not None) is refused
        if refused:
            assert "imaginary frequency of" in reason
            assert "should be zero" in reason
