# -*- coding: utf-8 -*-

"""One Hessian per species however many optimization blocks it takes
(Thermochemistry 5874, 2026-11-21): an optimization that does not converge in
'fragment max steps' continues in the same engine session, and only then are the
energy and the Hessian computed, at the final structure. Before, the second call
of optimize_and_hessian computed an analytic Hessian at the unconverged
structure and threw it away (~1 h for a 700-1000 basis-function target)."""

import contextlib

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from seamm_util import Q_  # noqa: E402

#: Methane, C-H 1.089 Å
CH4 = [
    ("C", (0.0, 0.0, 0.0)),
    ("H", (0.629, 0.629, 0.629)),
    ("H", (-0.629, -0.629, 0.629)),
    ("H", (-0.629, 0.629, -0.629)),
    ("H", (0.629, -0.629, -0.629)),
]
#: A harmonic potential E = 1/2 x.K.x (E_h, bohr) for the fake engine
K = np.diag(np.linspace(0.2, 1.6, 15))
BOHR = 0.529177210903  # Å
#: How far each optimization block moves every coordinate, Å
SHIFT = 0.01
MC = {
    "level": "ORCA:DFT@B3LYP/def2-SVP",
    "step": "ORCA",
    "method": "B3LYP",
    "options": {"mdi_capable": True, "analytic_hessian": True},
}


class Engine:
    """E = 1/2 x.K.x from where the structure is; counts its Hessians."""

    def __init__(self):
        self.x = np.zeros(15)
        self.n_hessians = 0

    def set_coordinates(self, x, units="Å"):
        self.x = np.asarray(x, dtype=float).reshape(-1) / BOHR

    def energy(self, units="kJ/mol"):
        e = 0.5 * self.x @ K @ self.x
        return e * 2625.4996394799 if units == "kJ/mol" else e

    def supports(self, command):
        return command == "<HESSIAN"

    def hessian(self):
        self.n_hessians += 1
        return K.copy()

    def fingerprints(self):
        return None


@pytest.fixture()
def methane():
    from molsystem import SystemDB

    db = SystemDB(filename="file:hessonce?mode=memory&cache=shared")
    configuration = db.create_system().create_configuration()
    configuration.atoms.append(
        x=[p[0] for _, p in CH4],
        y=[p[1] for _, p in CH4],
        z=[p[2] for _, p in CH4],
        symbol=[s for s, _ in CH4],
    )
    yield configuration
    db.close()


def _run(node, configuration, monkeypatch, blocks):
    """_compute_species_mc with the real driver over a fake warm engine whose
    optimization blocks converge as `blocks` says, each moving the structure."""
    import seamm_exec.hessian as H

    engine = Engine()
    script = list(blocks)
    sessions = []

    @contextlib.contextmanager
    def session(self, configuration, announce=True):
        xyz = np.array(configuration.atoms.get_coordinates(fractionals=False))
        s = H.EngineSession(self, engine, configuration, None, [6, 1, 1, 1, 1], xyz)

        def optimize(fmax=0.005, max_steps=200):
            converged = script.pop(0)
            s.coordinates = s.coordinates + SHIFT
            configuration.atoms.set_coordinates(s.coordinates, fractionals=False)
            self.last_optimization = {
                "converged": converged,
                "evaluations": 3,
                "max steps": max_steps,
            }
            return converged

        s.optimize = optimize
        sessions.append(s)
        yield s

    seen = {}

    def thermochemistry(_P, *args, **kwargs):
        """Stand-in: the energy from _E0, the thermal terms from where the
        vibrations' atoms are, so both show which structure they came from."""
        positions = node._vibrations.get_atoms().positions
        seen["positions"] = positions.copy()
        offset = float(np.mean(positions - np.array([p for _, p in CH4])))
        return {
            "total_energy": Q_(node._E0, "eV").m_as("E_h"),
            "ZPE": 100.0 + 1000.0 * offset,
            "H thermal": [[110.0 + 1000.0 * offset]],
        }

    monkeypatch.setattr(H.HessianDriver, "session", session)
    monkeypatch.setattr(node, "model_chemistry", lambda: MC)
    monkeypatch.setattr(node, "harmonic_thermochemistry", thermochemistry)
    node._hessian_route = "warm engine over MDI"
    node.model = MC["level"]
    said = []
    monkeypatch.setattr(node, "_say_continuing", lambda can, n: said.append(n))
    P = {
        "thermal corrections": "harmonic",
        "optimize fragments": "yes",
        "fragment convergence": Q_(0.005, "eV/Å"),
        "fragment max steps": 200,
        "step size": Q_(0.01, "Å"),
    }
    result = node._compute_species_mc(P, "C", configuration)
    assert script == [], "every block was run"
    return result, engine, sessions, said, seen


@pytest.mark.parametrize(
    "blocks, converged",
    [
        ([True], True),  # converged in the first block: as before
        ([False, True], True),  # continued once, then converged
        ([False, False], False),  # still not converged: flagged
    ],
)
def test_one_hessian_at_the_final_structure(
    node, methane, monkeypatch, blocks, converged
):
    result, engine, sessions, said, seen = _run(node, methane, monkeypatch, blocks)
    n = len(blocks)
    assert engine.n_hessians == 1
    assert len(sessions) == 1, "the continuation runs in the same engine session"
    assert said == [200] * (n - 1)
    assert result["converged"] is converged
    # The energy and the thermal terms are those of the final structure
    final = np.array([p for _, p in CH4]) + n * SHIFT
    assert np.allclose(seen["positions"], final)
    x = final.reshape(-1) / BOHR
    E_final = 0.5 * x @ K @ x * 2625.4996394799
    assert result["E"] == pytest.approx(E_final, rel=1e-9)
    assert result["ZPE"] == pytest.approx(100.0 + 1000.0 * n * SHIFT)
    assert result["H0"] == pytest.approx(result["E"] + result["ZPE"])
    assert result["H298"] == pytest.approx(result["E"] + 110.0 + 1000.0 * n * SHIFT)
    assert np.allclose(methane.atoms.get_coordinates(fractionals=False), final)
