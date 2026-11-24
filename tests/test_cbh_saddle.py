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


def _hessian(saddle, rotation=-0.05, k=40.0, curvature=-1.0):
    """A Cartesian Hessian (eV/Å^2) for water-shaped atoms: zero along the
    translations and rotations except for a spurious negative curvature
    `rotation` (eV/Å^2/amu) along one rotation, as real Hessians have (CaH2,
    20.1i), and `k` along the vibrations, with `curvature` along one vibration
    (mostly the oxygen moving along x) for a saddle point. Returns (Hessian,
    the saddle direction in Cartesian coordinates)."""
    from ase import Atoms

    atoms = Atoms(SYMBOLS, positions=START)
    m = np.repeat(atoms.get_masses(), 3)
    sqrt_m = np.sqrt(m)
    r = START - atoms.get_center_of_mass()
    q = []
    for axis in np.eye(3):
        q.append((np.tile(axis, 3) * sqrt_m))  # translations
    for axis in np.eye(3):
        q.append((np.cross(axis, r).ravel() * sqrt_m))  # rotations
    Q, _ = np.linalg.qr(np.array(q).T)
    Q = Q[:, :6]
    P = np.eye(9) - Q @ Q.T
    d = np.zeros(9)
    d[0] = 1.0  # the oxygen along x
    d = P @ d
    d /= np.linalg.norm(d)
    H = k * P
    if saddle:
        H += (curvature - k) * np.outer(d, d)
    rot = Q[:, 5] - Q[:, :3] @ (Q[:, :3].T @ Q[:, 5])
    rot /= np.linalg.norm(rot)
    H += rotation * np.outer(rot, rot)
    cartesian = sqrt_m[:, None] * H * sqrt_m[None, :]
    direction = (d / sqrt_m).reshape(3, 3)
    return cartesian, direction


def _vibrations(saddle, **kwargs):
    """VibrationsData for water-shaped atoms (see _hessian)."""
    from ase import Atoms
    from ase.vibrations import VibrationsData

    atoms = Atoms(SYMBOLS, positions=START)
    hessian, _ = _hessian(saddle, **kwargs)
    return VibrationsData(atoms, hessian.reshape(3, 3, 3, 3))


def _expected_frequency():
    """The saddle point's imaginary frequency (cm^-1), with no spurious
    rotation, from ASE."""
    frequencies = _vibrations(saddle=True, rotation=0.0).get_frequencies()
    return max(abs(f.imag) for f in frequencies)


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
    expected = _expected_frequency()
    assert expected > 2 * cbh_approach.IMAGINARY_TOLERANCE
    node._vibrations = _vibrations(saddle=True)
    frequency, mode = node._largest_imaginary_mode()
    assert frequency == pytest.approx(expected, rel=1e-6)
    # Along the saddle direction, the atom moving most moving 0.3 Å
    _, direction = _hessian(saddle=True)
    cosine = np.vdot(mode, direction) / np.linalg.norm(mode) / np.linalg.norm(direction)
    assert abs(cosine) == pytest.approx(1.0, abs=1e-6)
    assert np.max(np.linalg.norm(mode, axis=1)) == pytest.approx(
        cbh_approach.SADDLE_DISPLACEMENT
    )
    node._vibrations = _vibrations(saddle=False)
    assert node._largest_imaginary_mode() == (None, None)


def test_translations_and_rotations_are_not_saddle_modes(node):
    """A minimum whose near-zero rotation is imaginary above the tolerance
    (CaH2, ChemAI 5739: 20.1i among the six that should be zero) is not a saddle
    point."""
    from ase.units import invcm

    node._vibrations = _vibrations(saddle=False, rotation=-0.1)
    energies = node._vibrations.get_energies()
    assert max(abs(e.imag) for e in energies) / invcm > (
        cbh_approach.IMAGINARY_TOLERANCE
    )
    assert node._largest_imaginary_mode() == (None, None)


def test_zero_modes_of_linear_species_and_atoms():
    from ase import Atoms
    from ase.vibrations import VibrationsData

    def count(symbols, positions):
        atoms = Atoms(symbols, positions=positions)
        n = 3 * len(atoms)
        return cbh_approach._n_zero_modes(
            VibrationsData(atoms, np.eye(n).reshape(len(atoms), 3, len(atoms), 3))
        )

    assert count("Ca", [[0, 0, 0]]) == 3
    assert count("HCaH", [[0, 0, -2.0], [0, 0, 0], [0, 0, 2.0]]) == 5
    assert count(SYMBOLS, START) == 6


def _harmonic(node, configuration, minimum_when):
    """A stand-in for optimize + Hessian: a minimum when `minimum_when(x)` for
    the oxygen's x displacement; records the convergence asked for and the
    largest atomic displacement."""
    calls = []

    def harmonic(fmax):
        delta = configuration.atoms.xyz - START
        z = delta[0, 0]
        minimum = minimum_when(z)
        calls.append((fmax, z, np.max(np.linalg.norm(delta, axis=1))))
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
    assert calls[0][1] != 0.0
    assert all(
        moved == pytest.approx(cbh_approach.SADDLE_DISPLACEMENT)
        for _, _, moved in calls
    )
    tight = 0.005 * cbh_approach.SADDLE_TIGHTENING
    assert all(fmax == pytest.approx(tight) for fmax, _, _ in calls)
    assert result["E"] == -10.0 and result["imaginary"] == 0
    assert "saddle" not in result
    assert node._flags({"[CH3-]": result})["CBH saddle-point species"] == ""


def test_both_ways_are_tried_even_when_the_first_is_a_minimum(node):
    """Paul 2026-11-23: displace both ways and keep the lower minimum."""
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: True)
    result = node._escape_saddle(
        _P(), "[CH]=O", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    assert len(calls) == 2 and result["imaginary"] == 0
    record = result["saddle re-optimized"]
    assert record["minimum reached"] and len(record["rounds"]) == 1


def test_a_saddle_that_stays_is_flagged(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: False)
    result = node._escape_saddle(
        _P(), "F[C](F)F", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    assert len(calls) == 4  # two rounds, both ways
    assert len(result["saddle re-optimized"]["rounds"]) == 2
    assert not result["saddle re-optimized"]["minimum reached"]
    assert result["saddle"] == pytest.approx(_expected_frequency(), rel=1e-6)
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
        assert calls == [] and result["saddle"] == pytest.approx(
            _expected_frequency(), rel=1e-6
        )


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


# The displaced re-optimization as Paul set it (2026-11-23): both ways, the lower
# minimum kept, at most 'saddle rounds' rounds, recorded


class FakeProperties:
    def __init__(self):
        self.values = {}

    def exists(self, name):
        return name in self.values

    def add(self, name, *args, **kwargs):
        self.values.setdefault(name, None)

    def put(self, name, value):
        self.values[name] = value

    def get(self, pattern="*"):
        return {k: {"value": v} for k, v in self.values.items()}


def _configuration_with_properties():
    configuration = _configuration()
    configuration.properties = FakeProperties()
    return configuration


def test_the_lower_minimum_is_kept_and_restored(node):
    """Both directions reach a minimum; the first (+) is lower, so the
    configuration goes back to it: coordinates, Hessian and stored properties."""
    node.model = "Fake:M"
    configuration = _configuration_with_properties()
    node._vibrations = _vibrations(saddle=True)
    seen = []

    def harmonic(fmax):
        z = configuration.atoms.xyz[0, 0] - START[0, 0]
        E = -12.0 if not seen else -10.0  # the first attempt is lower
        seen.append(z)
        node._vibrations = _vibrations(saddle=False)
        configuration.properties.put("E#Thermochemistry#Fake:M", E)
        return {"E": E, "imaginary": 0}

    result = node._escape_saddle(
        _P(), "CC=C", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    assert len(seen) == 2 and result["E"] == -12.0
    assert configuration.atoms.xyz[0, 0] - START[0, 0] == pytest.approx(seen[0])
    assert configuration.properties.values["E#Thermochemistry#Fake:M"] == -12.0
    record = configuration.properties.values[
        "saddle re-optimization#Thermochemistry#Fake:M"
    ]
    assert record["rounds"][0]["chosen"] == "+"
    assert record["E chosen (kJ/mol)"] == -12.0 and record["minimum reached"]
    flags = node._flags({"CC=C": result, "C": {"E": 1.0}})
    assert flags["CBH saddle re-optimized species"] == "CC=C"
    assert flags["CBH saddle-point species"] == ""


def test_a_second_round_reaches_the_minimum(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    count = []

    def harmonic(fmax):
        count.append(fmax)
        minimum = len(count) == 3  # round 2, first direction
        node._vibrations = _vibrations(saddle=not minimum)
        return {"E": -float(len(count)), "imaginary": 0 if minimum else 1}

    result = node._escape_saddle(
        _P(), "CO", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    record = result["saddle re-optimized"]
    assert len(count) == 4 and len(record["rounds"]) == 2
    assert record["minimum reached"] and result["E"] == -3.0
    assert record["rounds"][0]["chosen"] == "-"  # the lower saddle (E -2)
    assert "saddle" not in result


@pytest.mark.parametrize(
    "changes, calls_made",
    [
        ({"saddle rounds": 1}, 2),
        ({"saddle threshold": Q_(10000.0, "1/cm")}, 0),  # below it: a minimum
    ],
)
def test_the_parameters(node, changes, calls_made):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: False)
    node._escape_saddle(
        _P(**changes), "CO", configuration, {"E": 0.0, "imaginary": 1}, harmonic
    )
    assert len(calls) == calls_made


def test_the_displacement_parameter(node):
    configuration = _configuration()
    node._vibrations = _vibrations(saddle=True)
    harmonic, calls = _harmonic(node, configuration, lambda z: True)
    node._escape_saddle(
        _P(**{"saddle displacement": Q_(0.1, "Å")}),
        "CO",
        configuration,
        {"E": 0.0, "imaginary": 1},
        harmonic,
    )
    assert all(moved == pytest.approx(0.1) for _, _, moved in calls)


# A real methyl rotation saddle: propene with its methyl rotated to the top of the
# barrier, with MMFF (RDKit) as the engine


def _mmff_engine(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    props = AllChem.MMFFGetMoleculeProperties(mol)
    return mol, props


def _mmff_hessian(mol, props, xyz, step=0.002):
    """The Cartesian Hessian (eV/Å^2) by central differences of MMFF energies.
    (RDKit's CalcGrad(positions) disagrees with its own energies by a factor of
    ~4.18 here, so the energies are used.)"""
    from rdkit.Chem import AllChem

    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
    flat = np.asarray(xyz, dtype=float).ravel()
    n3 = len(flat)

    def energy(*moves):
        x = flat.copy()
        for k, d in moves:
            x[k] += d
        return ff.CalcEnergy(x.tolist())

    E0 = energy()
    H = np.zeros((n3, n3))
    h = step
    for a in range(n3):
        H[a, a] = (energy((a, h)) - 2 * E0 + energy((a, -h))) / h**2
        for b in range(a):
            H[a, b] = H[b, a] = (
                energy((a, h), (b, h))
                - energy((a, h), (b, -h))
                - energy((a, -h), (b, h))
                + energy((a, -h), (b, -h))
            ) / (4 * h * h)
    return H * 0.0433641043  # kcal/mol/Å^2 -> eV/Å^2 (per molecule)


def _mmff_harmonic(node, configuration, mol, props, constraint=None):
    """optimize (MMFF, optionally with the methyl torsion held) + Hessian."""
    from ase import Atoms
    from ase.vibrations import VibrationsData
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    symbols = [a.GetSymbol() for a in mol.GetAtoms()]

    def harmonic(fmax):
        conf = mol.GetConformer()
        for i, (x, y, z) in enumerate(configuration.atoms.xyz):
            conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
        if constraint is not None:
            ff.MMFFAddTorsionConstraint(*constraint, False, 60.0, 60.0, 1.0e3)
        for _ in range(50):  # until converged (0)
            if ff.Minimize(maxIts=2000, forceTol=1e-8, energyTol=1e-14) == 0:
                break
        xyz = np.array(ff.Positions()).reshape(-1, 3)
        configuration.atoms.xyz = xyz
        E = AllChem.MMFFGetMoleculeForceField(mol, props).CalcEnergy(
            xyz.ravel().tolist()
        )
        H = _mmff_hessian(mol, props, xyz)
        n = len(symbols)
        node._vibrations = VibrationsData(
            Atoms(symbols, positions=xyz), H.reshape(n, 3, n, 3)
        )
        return {"E": E * 4.184, "imaginary": 0}

    return harmonic


def test_a_planted_methyl_saddle_ends_at_the_minimum(node):
    pytest.importorskip("ase")
    from rdkit.Chem import rdMolTransforms

    mol, props = _mmff_engine("CC=C")
    # An H on the methyl carbon (atom 0) and the torsion H-C0-C1=C2: the minimum
    # has a C-H eclipsing C=C (0, 120 or 240 deg); hold it at 60 deg, the top of the
    # barrier (MMFF: 8.44 against 6.48 kcal/mol)
    h = next(
        a.GetIdx() for a in mol.GetAtomWithIdx(0).GetNeighbors() if a.GetSymbol() == "H"
    )
    conf = mol.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, h, 0, 1, 2, 60.0)
    configuration = SimpleNamespace(n_atoms=mol.GetNumAtoms(), atoms=FakeAtoms())
    configuration.atoms.xyz = np.array(conf.GetPositions())
    configuration.atoms.symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    held = _mmff_harmonic(node, configuration, mol, props, constraint=(h, 0, 1, 2))
    saddle = held(0.001)  # optimized with the torsion held at the barrier top
    frequency, _ = node._largest_imaginary_mode(50.0)
    assert frequency is not None and frequency > 100.0  # the methyl torsion
    free = _mmff_harmonic(node, configuration, mol, props)
    result = node._escape_saddle(
        _P(**{"saddle threshold": Q_(50.0, "1/cm")}),
        "CC=C",
        configuration,
        dict(saddle),
        free,
    )
    record = result["saddle re-optimized"]
    assert record["minimum reached"]
    assert node._largest_imaginary_mode(50.0) == (None, None)
    barrier = saddle["E"] - result["E"]
    assert 3.0 < barrier < 15.0  # propene's methyl barrier, ~8 kJ/mol
    dihedral = rdMolTransforms.GetDihedralDeg(mol.GetConformer(), h, 0, 1, 2)
    # an H eclipses C=C again: the dihedral is a multiple of 120 degrees
    assert abs((dihedral + 60.0) % 120.0 - 60.0) < 10.0


def test_a_true_minimum_is_untouched(node):
    pytest.importorskip("ase")
    mol, props = _mmff_engine("CC=C")
    configuration = SimpleNamespace(n_atoms=mol.GetNumAtoms(), atoms=FakeAtoms())
    configuration.atoms.xyz = np.array(mol.GetConformer().GetPositions())
    configuration.atoms.symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    free = _mmff_harmonic(node, configuration, mol, props)
    minimum = free(0.001)
    xyz = configuration.atoms.xyz.copy()
    calls = []

    def counting(fmax):
        calls.append(fmax)
        return free(fmax)

    result = node._escape_saddle(_P(), "CC=C", configuration, dict(minimum), counting)
    assert calls == [] and result == minimum
    assert np.allclose(configuration.atoms.xyz, xyz)
