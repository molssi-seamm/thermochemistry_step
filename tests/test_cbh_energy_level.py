# -*- coding: utf-8 -*-

"""The CBH approach with the energies at another level than the geometries:
the energy level's single points (one batch through seamm_exec's Evaluator,
the target first), the combination E(energy level) + [H - E](geometry level),
the properties stored for the fragment store's harvest, and the reuse of
stored geometries and energies."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

import seamm_exec  # noqa: E402

from tests.conftest import ENERGIES  # noqa: E402
from tests.test_cbh_existing import _parameters  # noqa: E402

LOW = "Low:R2SCAN"
HIGH = {"level": "High:F12", "step": "High"}
# The geometry level's energies, offset from the high level's (kJ/mol)
OFFSET = {"CCF": 100.0, "C": 40.0, "CF": 70.0, "CC": 60.0}


class FakeEvaluator:
    """Stands in for seamm_exec.Evaluator: the 'energy level' energy of each
    configuration from ENERGIES, by its system's name."""

    calls = []

    def __init__(self, node, mc, **kwargs):
        assert mc is HIGH and kwargs["properties"] == ("energy",)
        self.submitted = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def submit(self, configuration, key):
        self.submitted.append((key, configuration))

    def results(self):
        from seamm_util.fingerprint import geometry_fingerprint

        names = []
        for key, configuration in self.submitted:
            name = configuration.system.name
            names.append(name)
            can = name[4:] if name.startswith("CBH:") else "CCF"
            yield seamm_exec.EvaluatorResult(
                key=key,
                ok=True,
                energy=ENERGIES[can],
                data={
                    "calculation": {
                        "program": "High",
                        "fingerprint": "fp:high",
                        "geometry fingerprint": geometry_fingerprint(
                            configuration.atoms.atomic_numbers,
                            configuration.atoms.get_coordinates(fractionals=False),
                        ),
                    }
                },
            )
        FakeEvaluator.calls.append(names)


@pytest.fixture()
def split(node, monkeypatch):
    """The node computing its species (stubbed: geometry-level energies with ZPE
    10 and H - E 15) and with the energy level HIGH."""
    FakeEvaluator.calls = []
    node.model = LOW
    monkeypatch.setattr(seamm_exec, "Evaluator", FakeEvaluator)
    monkeypatch.setattr(node, "_energy_level", lambda P: HIGH)
    computed = []

    def compute(P, can, generated, system=None, configuration=None, embed=False):
        computed.append(can)
        E = ENERGIES[can] + OFFSET[can]
        return {"E": E, "ZPE": 10.0, "H0": E + 10.0, "H298": E + 15.0}

    monkeypatch.setattr(node, "_compute_species", compute)
    node._computed = computed
    return node


def _computed(node, **changes):
    return _parameters(
        node,
        **{
            "cbh energies": "finite-difference harmonic",
            "thermal corrections": "harmonic",
            "energy model chemistry": "High:F12",
            **changes,
        },
    )


def test_energies_at_the_energy_level(split):
    data = split.run_cbh(_computed(split))
    # E(energy level) + [H - E](geometry level): the same reaction energy as the
    # high-level energies alone (CBH-1 conserves the molecules, so the constant
    # thermal part cancels), not the geometry level's
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert data["CBH energy level"] == "High:F12"
    assert data["CBH geometry level"] == LOW
    # The target first, alone; then the fragments in one batch
    assert FakeEvaluator.calls[0] == ["CBH:CCF"]
    assert sorted(FakeEvaluator.calls[1]) == ["CBH:C", "CBH:CC", "CBH:CF"]
    # Stored for the harvest
    configuration = split._find_fragment(split._test_db, "CF")
    p = configuration.properties
    assert (
        p.get("settings fingerprint#Thermochemistry#High:F12")[
            "settings fingerprint#Thermochemistry#High:F12"
        ]["value"]
        == "fp:high"
    )
    E = p.get("total energy#Thermochemistry#High:F12")[
        "total energy#Thermochemistry#High:F12"
    ]["value"]
    assert E == pytest.approx(-230.0 / 2625.4996394799)


def test_geometry_level_is_the_default(split, monkeypatch):
    """Without an energy level nothing runs on the Evaluator."""
    monkeypatch.setattr(split, "_energy_level", lambda P: None)
    data = split.run_cbh(
        _computed(split, **{"energy model chemistry": "the geometry level"})
    )
    assert FakeEvaluator.calls == []
    assert data["CBH energy level"] == LOW
    # The geometry level's own energies: the offsets change the reaction energy
    # by OFFSET[CF] + OFFSET[CC] - OFFSET[CCF] - OFFSET[C] = -10
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0 + 10.0, abs=1e-6)


def test_existing_results_mode_refuses_an_energy_level(split):
    with pytest.raises(RuntimeError, match="need the species computed"):
        split.run_cbh(
            _computed(split, **{"cbh energies": "existing results (run nothing)"})
        )


def test_stored_geometry_and_energy_are_reused(split, tmp_path):
    """A fragment whose geometry (at the geometry level) and energy (at the
    energy level, same settings fingerprint) are in the store is neither
    computed nor sent to the Evaluator."""
    from seamm_thermochemistry import EnergyRecord, FragmentStore, GeometryRecord

    path = tmp_path / "fragments.db"
    xyz = [
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.38],
        [1.03, 0.0, -0.36],
        [-0.51, 0.89, -0.36],
        [-0.51, -0.89, -0.36],
    ]
    with FragmentStore(path, create=True) as store:
        g = GeometryRecord.make(
            "CF",
            ["C", "F", "H", "H", "H"],
            xyz,
            E=ENERGIES["CF"] + OFFSET["CF"],
            ZPE=10.0,
            H=ENERGIES["CF"] + OFFSET["CF"] + 15.0,
            model_chemistry=LOW,
        )
        store.put_geometry(g)
        store.put_energy(
            EnergyRecord.make(
                g["geometry_hash"], ENERGIES["CF"], energy_fingerprint="fp:high"
            )
        )
    P = _computed(split, **{"fragment store": "use if available"})
    P["fragment store file"] = str(path)
    data = split.run_cbh(P)
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert "CF" not in split._computed
    assert all("CBH:CF" not in call for call in FakeEvaluator.calls)
    assert FakeEvaluator.calls[0] == ["CBH:CCF"]


def test_species_only(split, monkeypatch):
    """Only the target is computed and sent to the energy level, with the
    properties the harvest reads; no reactions and no ATcT (an ATcT that raises
    is never touched)."""
    import seamm_thermochemistry

    def no_atct(*args, **kwargs):
        raise AssertionError("ATcT must not be used")

    monkeypatch.setattr(seamm_thermochemistry, "ATcT", no_atct)
    P = _computed(
        split, **{"cbh calculation": "the species only (for the fragment store)"}
    )
    data = split.run_cbh(P)
    assert split._computed == ["CCF"]
    # The target alone, once
    assert len(FakeEvaluator.calls) == 1 and len(FakeEvaluator.calls[0]) == 1
    assert data == {
        "CBH energy level": "High:F12",
        "CBH geometry level": LOW,
        "CBH unconverged species": "",
        "CBH changed fragments": "",
    }
    assert "DfH 298K" not in data
    _, configuration = split.get_system_configuration()
    p = configuration.properties
    name = "settings fingerprint#Thermochemistry#High:F12"
    assert p.get(name)[name]["value"] == "fp:high"
    name = "total energy#Thermochemistry#High:F12"
    assert p.get(name)[name]["value"] == pytest.approx(
        ENERGIES["CCF"] / 2625.4996394799
    )


def test_species_only_refuses_existing_results(split):
    with pytest.raises(RuntimeError, match="computes the species"):
        split.run_cbh(
            _computed(
                split,
                **{
                    "cbh calculation": "the species only (for the fragment store)",
                    "cbh energies": "existing results (run nothing)",
                },
            )
        )


def test_species_only_description(split):
    from thermochemistry_step.cbh_approach import PRESETS

    text = split.describe_cbh(
        {
            "cbh calculation": "the species only (for the fragment store)",
            "preset": "default",
        }
    )
    assert "No CBH reactions" in text and PRESETS["default"]["energy"] in text
