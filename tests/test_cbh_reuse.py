# -*- coding: utf-8 -*-

"""Reuse of the CBH species (P5a of the fragment store, see seamm_thermochemistry's
docs/developer_guide/campaigns/2026-11-10):

* a species already computed in the job (an earlier iteration of a Loop) is
  taken from the job's database when its fingerprints are the expected ones;
* each species is looked up with the fingerprints the program predicts for IT,
  not the target's: ORCA gives a species with B or P (here: F, in the fake
  program) another fingerprint, so a B/P target's ordinary fragments must still
  find their stored energies;
* the target's fingerprints, not the last species', are the fallback;
* a computed fingerprint that is not the predicted one is reported.
"""

import types

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

import seamm_exec  # noqa: E402

from tests.conftest import ENERGIES  # noqa: E402
from tests.test_cbh_energy_level import (  # noqa: E402,F401
    HIGH,
    LOW,
    OFFSET,
    FakeEvaluator,
    _computed,
    split,
)

LOW_MC = {"level": LOW, "step": "Fake"}


class FakeProgram:
    """Predicts fingerprints as a program whose settings change with the
    elements: a species with F gets ':X' (as ORCA's NoCOSX for B and P)."""

    @staticmethod
    def expected_fingerprints(
        configuration, model_chemistry, *, kind, properties, options, version
    ):
        numbers = seamm_exec.structure_data(configuration)["atomic_numbers"]
        suffix = ":X" if 9 in numbers else ""
        return [f"fp:{model_chemistry['level']}{suffix}"]


def _fingerprint(configuration, level):
    numbers = seamm_exec.structure_data(configuration)["atomic_numbers"]
    return f"fp:{level}" + (":X" if 9 in numbers else "")


class ElementEvaluator(FakeEvaluator):
    """FakeEvaluator whose records carry the fake program's fingerprint for each
    species, and a version."""

    def results(self):
        for result in super().results():
            configuration = dict(self.submitted)[result.key]
            result.data["calculation"]["fingerprint"] = _fingerprint(
                configuration, HIGH["level"]
            )
            result.data["calculation"]["version"] = "1.0"
            yield result


@pytest.fixture()
def program(split, monkeypatch):  # noqa: F811
    """The split node with the fake program behind both levels."""
    split.flowchart = types.SimpleNamespace(
        plugin_manager=types.SimpleNamespace(get=lambda name: FakeProgram)
    )
    monkeypatch.setattr(split, "model_chemistry", lambda: LOW_MC)
    return split


def _xyz(configuration):
    return configuration.atoms.get_coordinates(fractionals=False)


def test_species_structure():
    from thermochemistry_step.cbh_approach import CBHApproach

    methyl = CBHApproach._species_structure("[CH3]")
    assert methyl.n_atoms == 4 and methyl.spin_multiplicity == 2
    hf = CBHApproach._species_structure("F")
    assert sorted(seamm_exec.structure_data(hf)["atomic_numbers"]) == [1, 9]


def test_a_targets_ordinary_fragments_use_their_own_fingerprint(
    program, monkeypatch, tmp_path
):
    """The target CCF has F, so its energy-level fingerprint is 'fp:High:F12:X';
    C and CC have none and are stored under 'fp:High:F12'. Looked up with the
    target's fingerprint they would miss and be recomputed (job 5684's boron and
    phosphorus targets); with their own they are found."""
    from seamm_thermochemistry import EnergyRecord, FragmentStore, GeometryRecord
    from seamm_thermochemistry.geometry_store import hash_geometry

    monkeypatch.setattr(seamm_exec, "Evaluator", ElementEvaluator)
    path = tmp_path / "fragments.db"
    db = program._test_db
    with FragmentStore(path, create=True) as store:
        for can in ("C", "CC", "CF"):
            configuration = program._find_fragment(db, can)
            symbols = configuration.atoms.symbols
            xyz = _xyz(configuration)
            E = ENERGIES[can] + OFFSET[can]
            g = GeometryRecord.make(
                can, symbols, xyz, E=E, ZPE=10.0, H=E + 15.0, model_chemistry=LOW
            )
            store.put_geometry(g)
            store.put_energy(
                EnergyRecord.make(
                    hash_geometry(symbols, xyz),
                    ENERGIES[can],
                    energy_fingerprint=_fingerprint(configuration, HIGH["level"]),
                )
            )
    FakeEvaluator.calls = []
    P = _computed(program, **{"fragment store": "use if available"})
    P["fragment store file"] = str(path)
    data = program.run_cbh(P)
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    # Only the target ran at the energy level: C, CC (no F) and CF (with F)
    # were all found under their own fingerprints
    assert FakeEvaluator.calls == [["CBH:CCF"]]
    assert program._computed == ["CCF"]


def test_without_a_prediction_the_targets_fingerprint_is_used(
    program, monkeypatch, tmp_path
):
    """A program that cannot predict: the target's fingerprint (today's rule),
    so the fragments without F miss and run."""
    from seamm_thermochemistry import EnergyRecord, FragmentStore, GeometryRecord
    from seamm_thermochemistry.geometry_store import hash_geometry

    program.flowchart = types.SimpleNamespace(
        plugin_manager=types.SimpleNamespace(get=lambda name: object())
    )
    monkeypatch.setattr(seamm_exec, "Evaluator", ElementEvaluator)
    path = tmp_path / "fragments.db"
    db = program._test_db
    with FragmentStore(path, create=True) as store:
        for can in ("C", "CC", "CF"):
            configuration = program._find_fragment(db, can)
            symbols, xyz = configuration.atoms.symbols, _xyz(configuration)
            E = ENERGIES[can] + OFFSET[can]
            store.put_geometry(
                GeometryRecord.make(
                    can, symbols, xyz, E=E, ZPE=10.0, H=E + 15.0, model_chemistry=LOW
                )
            )
            store.put_energy(
                EnergyRecord.make(
                    hash_geometry(symbols, xyz),
                    ENERGIES[can],
                    energy_fingerprint=_fingerprint(configuration, HIGH["level"]),
                )
            )
    FakeEvaluator.calls = []
    P = _computed(program, **{"fragment store": "use if available"})
    P["fragment store file"] = str(path)
    program.run_cbh(P)
    assert FakeEvaluator.calls[0] == ["CBH:CCF"]
    assert sorted(FakeEvaluator.calls[1]) == ["CBH:C", "CBH:CC"]


def _put(configuration, name, kind, value, units=None):
    p = configuration.properties
    if not p.exists(name):
        p.add(name, kind, units=units, description=name)
    p.put(name, value)


def _computed_in_job(configuration, can, frequency_fingerprint, frequencies=None):
    """The properties a computed species leaves on its configuration."""
    E = ENERGIES[can] + OFFSET[can]
    s = f"#Thermochemistry#{LOW}"
    _put(configuration, "total energy" + s, "float", E / 2625.4996394799, "E_h")
    _put(configuration, "ZPE" + s, "float", 10.0, "kJ/mol")
    _put(configuration, "H(298.15 K)" + s, "float", E + 15.0, "kJ/mol")
    _put(configuration, "frequency fingerprint" + s, "str", frequency_fingerprint)
    if frequencies is not None:
        _put(configuration, "frequencies" + s, "json", frequencies)
    s = f"#Thermochemistry#{HIGH['level']}"
    _put(
        configuration,
        "total energy" + s,
        "float",
        ENERGIES[can] / 2625.4996394799,
        "E_h",
    )
    _put(configuration, "settings fingerprint" + s, "str", "fp:high")
    # The record names the geometry the single point ran on: this one
    from seamm_util.fingerprint import geometry_fingerprint

    ghash = geometry_fingerprint(
        configuration.atoms.atomic_numbers, _xyz(configuration)
    )
    _put(
        configuration,
        "calculation" + s,
        "json",
        {"fingerprint": "fp:high", "geometry fingerprint": ghash},
    )


@pytest.fixture()
def target_sets_fallback(split, monkeypatch):  # noqa: F811
    """Computing the target sets the fallback fingerprints, as
    _store_fingerprints does (no prediction: the fake program has none)."""
    compute = split._compute_species

    def computing(P, can, *args, **kwargs):
        if can == "CCF":
            split._settings_fingerprint = "fp:low"
            split._frequency_fingerprint = "fp:low"
        return compute(P, can, *args, **kwargs)

    monkeypatch.setattr(split, "_compute_species", computing)
    return split


@pytest.mark.parametrize(
    "fingerprint, frequencies, reused",
    [
        ("fp:low", [300.0, 1500.0], True),
        ("fp:other", [300.0, 1500.0], False),  # other settings: computed
        ("fp:low", [-50.0, 300.0], False),  # an imaginary mode: computed
        ("fp:low", [-15.0, 300.0], True),  # below 20i: noise
    ],
)
def test_species_already_in_the_job(
    target_sets_fallback, fingerprint, frequencies, reused
):
    node = target_sets_fallback
    configuration = node._find_fragment(node._test_db, "C")
    _computed_in_job(configuration, "C", fingerprint, frequencies)
    FakeEvaluator.calls = []
    data = node.run_cbh(_computed(node))
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    sent = [name for call in FakeEvaluator.calls for name in call]
    if reused:
        # Neither the geometry level nor the energy level ran for methane
        assert "C" not in node._computed and "CBH:C" not in sent
    else:
        # The geometry level runs again. (The stub leaves the coordinates as they
        # were, so the energy level's single point on them is still valid; a
        # moved geometry is tested below.)
        assert "C" in node._computed


def test_the_fallback_is_the_targets_fingerprint(node):
    """The first species (the target) sets the fallback, the version and the
    Hessian route; a later species (a boron fragment, say) does not replace
    them, or every later lookup would use its fingerprint."""
    node.model = "Test:PM7"
    node._reset_fingerprints()
    configuration = node._test_db.system.configuration
    first = {"fingerprint": "fp:target", "version": "6.1.1"}
    node._fingerprints = {"energy": first, "hessian": None}
    node._store_fingerprints(configuration, "finite differences over MDI", "CCF")
    node._fingerprints = {
        "energy": {"fingerprint": "fp:boron", "version": "6.1.1"},
        "hessian": None,
    }
    node._store_fingerprints(configuration, "analytic", "B")
    assert node._settings_fingerprint == "fp:target"
    assert node._frequency_fingerprint == "fp:target"
    assert node._fingerprint_version == "6.1.1"
    assert node._hessian_how == "finite differences over MDI"
    node._reset_fingerprints()
    assert node._settings_fingerprint is None and node._fingerprint_version is None


def test_the_self_check_reports_a_drift(program, caplog):
    """A computed fingerprint that is not the predicted one is reported."""
    program._reset_fingerprints()
    program.model = LOW
    configuration = program._find_fragment(program._test_db, "CF")
    program._fingerprints = {
        "energy": {"fingerprint": f"fp:{LOW}:X", "version": "1.0"},
        "hessian": None,
    }
    program._store_fingerprints(configuration, None, "CF")
    assert "not the one predicted" not in caplog.text
    program._fingerprints = {
        "energy": {"fingerprint": "fp:drifted", "version": "1.0"},
        "hessian": None,
    }
    program._store_fingerprints(configuration, None, "CF")
    assert "not the one predicted" in caplog.text


def test_an_energy_on_another_geometry_is_not_reused(target_sets_fallback):
    """The species is reused at the geometry level, but its energy-level record
    names another geometry (it was computed again in place): the single point
    runs again."""
    node = target_sets_fallback
    configuration = node._find_fragment(node._test_db, "C")
    _computed_in_job(configuration, "C", "fp:low", [300.0])
    name = f"calculation#Thermochemistry#{HIGH['level']}"
    configuration.properties.put(
        name, {"fingerprint": "fp:high", "geometry fingerprint": "sha256:moved"}
    )
    FakeEvaluator.calls = []
    node.run_cbh(_computed(node))
    sent = [n for call in FakeEvaluator.calls for n in call]
    assert "C" not in node._computed and "CBH:C" in sent
