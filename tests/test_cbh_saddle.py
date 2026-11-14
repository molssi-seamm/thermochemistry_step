# -*- coding: utf-8 -*-

"""Species that optimize to a saddle point are displaced along the imaginary
mode and optimized again (P5 item 6: 15 saddle fragments in PPE1694 phase 1,
e.g. [CH]=O 1528i), and species whose starting structure is broken are not
computed (bridged Al2Cl2 from UFF, PPE1694 5729)."""

from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from seamm_util import Q_  # noqa: E402

from thermochemistry_step import cbh_approach  # noqa: E402

SYMBOLS = "OHH"
START = np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]])


def _vibrations(saddle):
    """VibrationsData for water-shaped atoms whose Hessian is diagonal, with a
    negative curvature along the oxygen's z for a saddle point."""
    from ase import Atoms
    from ase.vibrations import VibrationsData

    atoms = Atoms(SYMBOLS, positions=START)
    diagonal = np.ones(9)
    if saddle:
        diagonal[2] = -0.5  # eV/Å^2: about 92i cm^-1
    return VibrationsData(atoms, np.diag(diagonal).reshape(3, 3, 3, 3))


class FakeAtoms:
    def __init__(self):
        self.xyz = START.copy()
        self.symbols = list(SYMBOLS)

    def get_coordinates(self, fractionals=True):
        return self.xyz.tolist()

    def set_coordinates(self, xyz, fractionals=True):
        self.xyz = np.array(xyz, dtype=float)


def _configuration():
    return SimpleNamespace(n_atoms=3, atoms=FakeAtoms())


def _P(**changes):
    return {
        "saddle points": "displace and re-optimize",
        "optimize fragments": "yes",
        "fragment convergence": Q_(0.005, "eV/Å"),
        **changes,
    }


def test_the_largest_imaginary_mode(node):
    node._vibrations = _vibrations(saddle=True)
    frequency, mode = node._largest_imaginary_mode()
    assert frequency == pytest.approx(92, abs=3)
    # Along the oxygen's z only, its largest atomic displacement 0.3 Å
    assert np.abs(mode[0, 2]) == pytest.approx(cbh_approach.SADDLE_DISPLACEMENT)
    assert np.abs(mode).sum() == pytest.approx(cbh_approach.SADDLE_DISPLACEMENT)
    node._vibrations = _vibrations(saddle=False)
    assert node._largest_imaginary_mode() == (None, None)


def _harmonic(node, configuration, minimum_when):
    """A stand-in for optimize + Hessian: a minimum when `minimum_when(z)` for
    the oxygen's z after displacement; records the convergence asked for."""
    calls = []

    def harmonic(fmax):
        z = configuration.atoms.xyz[0, 2]
        minimum = minimum_when(z)
        calls.append((fmax, z))
        node._vibrations = _vibrations(saddle=not minimum)
        return {"E": -10.0 if minimum else 0.0, "imaginary": 0 if minimum else 1}

    return harmonic, calls


def test_a_saddle_is_displaced_to_a_minimum(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    saddle = {"E": 0.0, "imaginary": 1}
    # The minimum lies on the side opposite the first displacement (whichever
    # sign the diagonalization gave the mode)
    first = []

    def other_side(z):
        first.append(z)
        return z * first[0] < 0

    harmonic, calls = _harmonic(node, configuration, other_side)
    result = node._escape_saddle(_P(), "[CH3-]", configuration, saddle, harmonic)
    assert len(calls) == 2  # + stayed a saddle, - reached the minimum
    assert calls[0][1] == pytest.approx(-calls[1][1])
    assert abs(calls[0][1]) == pytest.approx(cbh_approach.SADDLE_DISPLACEMENT)
    tight = 0.005 * cbh_approach.SADDLE_TIGHTENING
    assert all(fmax == pytest.approx(tight) for fmax, _ in calls)
    assert result["E"] == -10.0 and result["imaginary"] == 0
    assert "saddle" not in result
    assert node._flags({"[CH3-]": result})["CBH saddle-point species"] == ""


def test_one_try_is_enough_when_it_reaches_a_minimum(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: True)
    result = node._escape_saddle(
        _P(), "[CH]=O", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    assert len(calls) == 1 and result["imaginary"] == 0


def test_a_saddle_that_stays_is_flagged(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: False)
    result = node._escape_saddle(
        _P(), "F[C](F)F", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    assert len(calls) == 2
    assert result["saddle"] == pytest.approx(92, abs=3)
    flags = node._flags({"F[C](F)F": result, "C": {"E": 1.0}})
    assert flags["CBH saddle-point species"] == "F[C](F)F"


def test_keep_and_flag_does_not_displace(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: True)
    saddle = {"E": 0.0, "imaginary": 1}
    for P in (
        _P(**{"saddle points": "keep and flag"}),
        _P(**{"optimize fragments": "no"}),
    ):
        result = node._escape_saddle(P, "[CH2]O", configuration, dict(saddle), harmonic)
        assert calls == [] and result["saddle"] == pytest.approx(92, abs=3)


def test_a_minimum_is_left_alone(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=False)
    harmonic, calls = _harmonic(node, configuration, lambda z: True)
    result = node._escape_saddle(_P(), "C", configuration, {"E": 1.0}, harmonic)
    assert calls == [] and result == {"E": 1.0}


def _molsystem_configuration(tmp_path, smiles):
    from molsystem import SystemDB

    db = SystemDB(filename=f"file:{tmp_path / 'start.db'}")
    system = db.create_system(name=f"CBH:{smiles}")
    configuration = system.create_configuration(name="initial")
    configuration.from_smiles(smiles)
    return configuration


@pytest.mark.parametrize(
    "smiles", ["[AlH2-]1[Cl+][AlH2-][Cl+]1", "FS(F)(F)(F)(F)F", "F[P-](F)(F)(F)(F)F"]
)
def test_force_fields_do_not_rebuild_metal_or_hypervalent_species(
    tmp_path, node, smiles
):
    """The 'lowest MMFF conformer' start keeps molsystem's structure for species
    the force fields break; it is sane, so it is computed."""
    configuration = _molsystem_configuration(tmp_path, smiles)
    before = np.array(configuration.atoms.get_coordinates(fractionals=False))
    node._embed_lowest(configuration, smiles, n=5)
    after = np.array(configuration.atoms.get_coordinates(fractionals=False))
    assert np.allclose(before, after)
    node._check_start(smiles, configuration)  # does not raise


def test_a_broken_start_is_not_computed(tmp_path, node):
    configuration = _molsystem_configuration(tmp_path, "[AlH2-]1[Cl+][AlH2-][Cl+]1")
    xyz = np.array(configuration.atoms.get_coordinates(fractionals=False))
    chlorines = [i for i, s in enumerate(configuration.atoms.symbols) if s == "Cl"]
    xyz[chlorines[1]] = xyz[chlorines[0]]
    configuration.atoms.set_coordinates(xyz, fractionals=False)
    with pytest.raises(RuntimeError, match="starting structure .* is broken"):
        node._check_start("[AlH2-]1[Cl+][AlH2-][Cl+]1", configuration)
