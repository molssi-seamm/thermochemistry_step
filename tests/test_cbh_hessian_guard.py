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
    """Analytic Hessians as listed; forces (for finite differences) from the
    potential 1/2 x.F.x."""

    def __init__(self, hessians, F=GOOD):
        self.x = np.zeros(15)
        self.hessians = list(hessians)
        self.F = F
        self.n_hessians = 0

    def forces(self, units="hartree/bohr"):
        return -(self.F @ self.x)

    def set_coordinates(self, x, units="Å"):
        x = np.asarray(x, dtype=float).reshape(-1)
        self.x = x if units == "bohr" else x / BOHR

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


def _run(node, configuration, monkeypatch, hessians, F=GOOD):
    import seamm_exec.hessian as H

    engine = Engine(hessians, F)
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
    result, engine, sessions, _ = _run(node, methane("C"), monkeypatch, [BAD])
    # once more, at the same structure, in a new engine session, by finite
    # differences of the forces (the analytic one is not asked for again)
    assert engine.n_hessians == 1 and len(sessions) == 2
    assert "hessian refused" not in result and result["ZPE"] == 100.0
    assert np.allclose(
        node._vibrations.get_hessian_2d(),
        Q_(GOOD, "hartree/bohr**2").m_as("eV/Å**2"),
        atol=1e-4,
    )


def test_a_target_whose_hessian_stays_broken_fails(node, methane, monkeypatch):
    with pytest.raises(RuntimeError, match="Hessian of C is broken"):
        _run(node, methane("C"), monkeypatch, [BAD], F=BAD)


def test_a_fragment_whose_hessian_stays_broken_is_flagged(node, methane, monkeypatch):
    result, engine, sessions, displaced = _run(
        node, methane("CBH:C"), monkeypatch, [BAD], F=BAD
    )
    assert engine.n_hessians == 1 and len(sessions) == 2
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


def test_a_replacement_with_a_point_in_another_state_is_refused(
    node, methane, monkeypatch
):
    """The replacement finite-difference Hessian is made with every point checked
    for a state change and tested against the gradients (seamm-exec 2026.11.26);
    a failure refuses the species (science 2026-11-26)."""
    import seamm_exec.hessian as H

    monkeypatch.setattr(
        H, "check_fd_points", lambda *a, **k: ["coordinate 4 +0.01 bohr: <S**2> moved"]
    )
    result, engine, sessions, displaced = _run(
        node, methane("CBH:C"), monkeypatch, [BAD]
    )
    assert "S**2" in result["hessian refused"]
    assert "ZPE" not in result and displaced == []


def test_a_sound_replacement_is_checked_against_the_gradients(
    node, methane, monkeypatch
):
    result, engine, sessions, _ = _run(node, methane("C"), monkeypatch, [BAD])
    assert "hessian refused" not in result
    check = node._hessian_check
    assert check is not None and not check["refused"] and len(check["errors"]) == 3
    assert node._fd_info["tight SCF"] is False  # this engine cannot tighten
    assert len(node._fd_info["zero modes"]) == 6
