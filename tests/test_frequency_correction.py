# -*- coding: utf-8 -*-

"""Per-frequency corrections (PM7: Rozanska et al. 2014, eq 1) in the harmonic
thermochemistry, the ZPE scale factor and the fingerprints."""

from ase import Atoms, units
from ase.vibrations import VibrationsData
import numpy as np
import pytest

from seamm_util import Q_

#: PM7 CH4 frequencies (cm^-1), from the step's own run (SEAMM_DEV, 2026-11-13)
CH4_PM7 = [1290.3, 1290.4, 1306.3, 1306.5, 1306.9, 2688.3, 2688.6, 2688.8, 2814.5]
CH4_XYZ = [
    [0.0, 0.0, 0.0],
    [0.629, 0.629, 0.629],
    [-0.629, -0.629, 0.629],
    [-0.629, 0.629, -0.629],
    [0.629, -0.629, -0.629],
]


def _vibrations(symbols, positions, frequencies):
    """An ASE VibrationsData whose real frequencies are exactly `frequencies`
    (cm^-1): a Hessian built in the space orthogonal to the translations and
    rotations."""
    atoms = Atoms(symbols, positions=positions)
    n = len(atoms)
    masses = atoms.get_masses()
    r = atoms.positions - atoms.get_center_of_mass()
    vectors = []
    for k in range(3):
        v = np.zeros((n, 3))
        v[:, k] = 1.0
        vectors.append((v * np.sqrt(masses)[:, None]).ravel())
        v = np.cross(np.eye(3)[k], r)
        vectors.append((v * np.sqrt(masses)[:, None]).ravel())
    U, S, _ = np.linalg.svd(np.array(vectors).T, full_matrices=True)
    rank = int(np.sum(S > 1e-6))
    null = U[:, rank:]
    assert null.shape[1] == len(frequencies)
    scale = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
    eigenvalues = (np.asarray(frequencies) * units.invcm / scale) ** 2
    hessian_mw = null @ np.diag(eigenvalues) @ null.T
    sqm = np.repeat(np.sqrt(masses), 3)
    hessian = hessian_mw * np.outer(sqm, sqm)
    return atoms, VibrationsData(atoms, hessian.reshape(n, 3, n, 3))


def _P(**extra):
    P = {
        "# imaginary modes": "perceive",
        "print frequencies": "all",
        "spin multiplicity": "from system",
        "symmetry number": "from system",
        "frequency correction": "from the table",
        "frequency correction range": "full",
    }
    P.update(extra)
    return P


@pytest.fixture()
def harmonic(node, monkeypatch):
    """The node set up for a harmonic analysis of a structure."""
    from molsystem import SystemDB

    db = SystemDB(filename="file:fc?mode=memory&cache=shared")
    cited = []
    monkeypatch.setattr(
        type(node), "_cite_frequency_correction", lambda self, c: cited.append(c.id)
    )

    def setup(symbols, positions, frequencies, model="MOPAC:SQM@PM7"):
        configuration = db.create_system().create_configuration()
        configuration.atoms.append(
            x=[p[0] for p in positions],
            y=[p[1] for p in positions],
            z=[p[2] for p in positions],
            symbol=list(symbols),
        )
        _, vibrations = _vibrations(symbols, positions, frequencies)
        node._working_configuration = configuration
        node._vibrations = vibrations
        node._E0 = 0.0
        node.model = model
        return node

    setup.cited = cited
    yield setup
    db.close()


def _zpe(frequencies):
    return Q_(0.5 * sum(frequencies), "1/cm").m_as("kJ/mol")


def test_ch4_zpe_with_and_without_the_correction(harmonic):
    from seamm_thermochemistry import frequency_correction

    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    corrected = frequency_correction("MOPAC:SQM@PM7")(CH4_PM7)

    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    assert data["vibrational frequencies"] == pytest.approx(CH4_PM7, abs=0.01)
    assert data["corrected vibrational frequencies"] == pytest.approx(
        corrected, abs=0.02
    )
    assert data["ZPE"] == pytest.approx(_zpe(corrected), abs=0.01)
    assert data["ZPE uncorrected"] == pytest.approx(_zpe(CH4_PM7), abs=0.01)
    assert data["ZPE"] - data["ZPE uncorrected"] == pytest.approx(6.19, abs=0.02)
    assert data["frequency correction"] == "fc:PM7-Rozanska2014:v1:full"
    assert data["frequency correction status"] == "applied"
    assert harmonic.cited == ["PM7-Rozanska2014"]

    plain = node.harmonic_thermochemistry(
        _P(**{"frequency correction": "none"}), [298.15], "K", [1.0], "bar"
    )
    assert plain["ZPE"] == pytest.approx(_zpe(CH4_PM7), abs=0.01)
    assert "frequency correction" not in plain
    # The entropy comes from the corrected frequencies too: the bends drop by
    # ~20 cm^-1 and dominate it, so S rises slightly
    assert data["S"][0][0] - plain["S"][0][0] == pytest.approx(0.049, abs=0.005)


def test_high_frequency_only(harmonic):
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    P = _P(**{"frequency correction range": "high-frequency only"})
    data = node.harmonic_thermochemistry(P, [298.15], "K", [1.0], "bar")
    assert data["corrected vibrational frequencies"][:5] == pytest.approx(
        CH4_PM7[:5], abs=0.01
    )
    assert data["frequency correction"].endswith(":high-frequency only")


def test_a_mode_above_the_fit_is_left_alone(harmonic):
    """PM7's H2 stretch, 5577 cm^-1, would become 14402 cm^-1."""
    node = harmonic(["H", "H"], [[0, 0, 0], [0, 0, 0.74]], [5577.0])
    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    assert data["corrected vibrational frequencies"] == pytest.approx([5577.0], abs=0.1)
    assert data["frequency correction status"] == "partly not applicable"
    assert data["ZPE"] == pytest.approx(data["ZPE uncorrected"])


def test_no_correction_for_other_mopac_methods(harmonic):
    node = harmonic(
        ["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7, model="MOPAC:SQM@PM6-ORG"
    )
    correction, _, note = node.frequency_correction(_P())
    assert correction is None
    assert "No frequency correction is available for MOPAC:SQM@PM6-ORG" in note
    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    assert "frequency correction" not in data
    assert data["ZPE"] == pytest.approx(_zpe(CH4_PM7), abs=0.01)


def test_the_discontinuity_is_reported(harmonic):
    """A raw frequency next to the 2500 cm^-1 switch is flagged."""
    from thermochemistry_step.thermochemistry import Thermochemistry  # noqa: F401

    freqs = [1290.3, 1290.4, 1306.3, 1306.5, 1306.9, 2490.0, 2688.6, 2688.8, 2814.5]
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, freqs)
    from seamm_thermochemistry import frequency_correction

    correction = frequency_correction("MOPAC:SQM@PM7")
    atoms = node.ase_atoms(node._working_configuration)
    text = node._frequency_correction_text(
        _P(), correction, "full", freqs, correction(freqs), atoms, 0
    )
    assert "discontinuous or not validated" in text and "2490.0" in text


class _Modes:
    """Stands in for VibrationsData: the given modes, in order."""

    def __init__(self, modes):
        self.modes = np.array(modes)

    def get_energies_and_modes(self):
        return np.ones(len(self.modes)), self.modes


def test_xh_stretches_outside_the_fit_set():
    from thermochemistry_step.thermochemistry import _xh_stretches

    r = 1.48
    xyz = [[0, 0, 0], [r, 0, 0], [-r, 0, 0], [0, r, 0], [0, -r, 0]]
    atoms = Atoms("SiH4", positions=xyz)
    stretch = np.zeros((5, 3))
    stretch[1] = [1.0, 0, 0]  # H1 along its Si-H bond
    bend = np.zeros((5, 3))
    bend[1] = [0, 1.0, 0]  # H1 across it
    symmetric = np.zeros((5, 3))
    for i, (x, y, _) in enumerate(xyz[1:], start=1):
        symmetric[i] = [x / r, y / r, 0]
    vibrations = _Modes([stretch, bend, symmetric])
    found = _xh_stretches(
        vibrations, 0, [2310.0, 900.0, 2320.0], atoms, ("C", "H", "N", "O")
    )
    assert found == [(2310.0, "Si-H"), (2320.0, "Si-H")]
    # Carbon is in the fit set: nothing to report for CH4
    ch4 = Atoms(
        "CH4",
        positions=[[0, 0, 0], [1.09, 0, 0], [-1.09, 0, 0], [0, 1.09, 0], [0, -1.09, 0]],
    )
    assert (
        _xh_stretches(vibrations, 0, [2310.0, 900.0, 2320.0], ch4, ("C", "H", "N", "O"))
        == []
    )


def test_the_correction_and_a_zpe_factor_exclude_each_other(node):
    node.model = "MOPAC:SQM@PM7"
    node._set_frequency_correction(_P(**{"thermal corrections": "harmonic"}))
    assert node._frequency_correction_tag == "fc:PM7-Rozanska2014:v1:full"
    scale, source, _ = node._zpe_scale({"ZPE scale factor": "from the table"})
    assert scale == 1.0 and "corrected" in source
    with pytest.raises(ValueError, match="cannot also be scaled"):
        node._zpe_scale({"ZPE scale factor": "0.98"})


def test_inert_without_thermal_corrections(node):
    node.model = "MOPAC:SQM@PM7"
    node._set_frequency_correction(
        _P(**{"thermal corrections": "none (the energy is the enthalpy, e.g. MOPAC)"})
    )
    assert node._frequency_correction_tag is None
    assert node._corrected_fingerprint("abc") == "abc"


def test_the_frequency_fingerprint_carries_the_correction(node):
    node.model = "MOPAC:SQM@PM7"
    node._set_frequency_correction(_P(**{"thermal corrections": "harmonic"}))
    assert node._corrected_fingerprint("abc") == "abc+fc:PM7-Rozanska2014:v1:full"
    assert node._corrected_fingerprint("") == ""
    node._set_frequency_correction(_P(**{"frequency correction": "none"}))
    assert node._corrected_fingerprint("abc") == "abc"
