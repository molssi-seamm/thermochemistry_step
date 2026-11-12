# -*- coding: utf-8 -*-

"""A composite recipe that combines several calculations as the CBH energy level:
the components run through composite_step's CompositeEvaluator (each as its own
program's batch), are looked up per component in the job and the fragment store,
never as a composite, and are stored as components; the composite record goes
with the results under its own property names."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")
pytest.importorskip("composite_step")

import seamm_exec  # noqa: E402
from composite_step import evaluate as E  # noqa: E402
from composite_step import recipes as R  # noqa: E402

from tests.conftest import ENERGIES  # noqa: E402
from tests.test_cbh_energy_level import (  # noqa: E402,F401
    LOW,
    OFFSET,
    _computed,
    split,
)

A, B = "High:QC@A/x", "High:QC@B/y"
SHIFT = 7.0
RECIPE = R.Recipe.from_dict(
    {
        "recipe": "toy",
        "version": 1,
        "geometry": {"level": "Low:DFT@R2SCAN"},
        "components": {"a": {"level": A}, "b": {"level": B}},
        "combine": "E0 = a - b\n",
    }
)
MC = {"level": "SEAMM:CMP@toy", "method": "toy", "step": "Composite", "options": {}}


class ComponentEvaluator:
    """Stands in for seamm_exec.Evaluator for each component: a = E + SHIFT and
    b = SHIFT, so the composite E0 = a - b is the conftest's ENERGIES."""

    calls = []
    fail = None

    def __init__(self, node, mc, **kwargs):
        assert kwargs["properties"] == ("energy",)
        self.mc = mc
        self.submitted = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def submit(self, configuration, key, options=None):
        self.submitted.append((key, configuration))

    def results(self):
        names = []
        for key, configuration in self.submitted:
            name = configuration.system.name
            names.append(name)
            can = name[4:] if name.startswith("CBH:") else "CCF"
            method = self.mc["method"]
            if ComponentEvaluator.fail == (method, can):
                yield seamm_exec.EvaluatorResult(
                    key=key, ok=False, reason="SCF did not converge"
                )
                continue
            energy = ENERGIES[can] + SHIFT if method == "A/x" else SHIFT
            yield seamm_exec.EvaluatorResult(
                key=key,
                ok=True,
                energy=energy,
                data={
                    "calculation": {"program": "High", "fingerprint": f"fp:{method}"}
                },
            )
        ComponentEvaluator.calls.append((self.mc["level"], names))


@pytest.fixture()
def composite(split, monkeypatch):  # noqa: F811
    ComponentEvaluator.calls = []
    ComponentEvaluator.fail = None
    monkeypatch.setattr(seamm_exec, "Evaluator", ComponentEvaluator)

    def resolve(level, available=None):
        method = level.split("@")[1]
        return {"level": level, "method": method, "step": "High", "options": {}}

    monkeypatch.setattr(E, "resolve_level", resolve)

    def energy_level(P):
        split._energy_recipe = RECIPE
        return MC

    monkeypatch.setattr(split, "_energy_level", energy_level)
    monkeypatch.setattr(
        split,
        "_expected_fingerprints",
        lambda structure, mc, kind, version=None: (
            None if mc is None else [f"fp:{mc['method']}"]
        ),
    )
    return split


def _property(configuration, name):
    return configuration.properties.get(name)[name]["value"]


def test_composite_energy_level(composite):
    data = composite.run_cbh(_computed(composite))
    # E0 = a - b = the conftest's energies: the same reaction as one energy level
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert data["CBH energy level"] == "SEAMM:CMP@toy"
    # Each component ran every species once, as its own batch
    levels = sorted(level for level, _ in ComponentEvaluator.calls)
    assert levels == [A, B]
    for _, names in ComponentEvaluator.calls:
        assert sorted(names) == ["CBH:C", "CBH:CC", "CBH:CCF", "CBH:CF"]
    # Components stored for the harvest under their own levels; the composite
    # only under "composite ..." names, which the harvest does not read
    cf = composite._find_fragment(composite._test_db, "CF")
    assert _property(cf, f"total energy#Thermochemistry#{A}") == pytest.approx(
        (ENERGIES["CF"] + SHIFT) / 2625.4996394799
    )
    assert _property(cf, f"settings fingerprint#Thermochemistry#{B}") == "fp:B/y"
    assert _property(cf, "composite energy#Thermochemistry#SEAMM:CMP@toy") == (
        pytest.approx(ENERGIES["CF"] / 2625.4996394799)
    )
    record = _property(cf, "composite calculation#Thermochemistry#SEAMM:CMP@toy")
    assert record["components"]["a"]["fingerprint"] == "fp:A/x"
    assert not cf.properties.exists("total energy#Thermochemistry#SEAMM:CMP@toy")
    assert not cf.properties.exists(
        "settings fingerprint#Thermochemistry#SEAMM:CMP@toy"
    )


def test_a_failed_component_names_itself(composite):
    ComponentEvaluator.fail = ("B/y", "CF")
    with pytest.raises(RuntimeError, match=r"CF failed: the b component"):
        composite.run_cbh(_computed(composite))


def test_store_lookup_is_by_component_never_composite(composite, tmp_path, monkeypatch):
    """A stored component row is reused for that component alone; the store is
    never asked for the composite."""
    from seamm_thermochemistry import (
        EnergyRecord,
        FragmentStore,
        GeometryRecord,
    )

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
                g["geometry_hash"], ENERGIES["CF"] + SHIFT, energy_fingerprint="fp:A/x"
            )
        )
    asked = []
    original = FragmentStore.energy

    def energy(self, *args, **kwargs):
        asked.append(kwargs.get("model_chemistry"))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(FragmentStore, "energy", energy)
    P = _computed(composite, **{"fragment store": "use if available"})
    P["fragment store file"] = str(path)
    data = composite.run_cbh(P)
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert asked and not any(str(m).startswith("SEAMM:CMP@") for m in asked)
    # CF's component a came from the store, with its fingerprint in the
    # composite record; its component b was computed
    cf = composite._find_fragment(composite._test_db, "CF")
    record = _property(cf, "composite calculation#Thermochemistry#SEAMM:CMP@toy")
    assert record["components"]["a"]["fingerprint"] == "fp:A/x"
    calls = dict(ComponentEvaluator.calls)
    assert "CBH:CF" not in calls[A]
    assert "CBH:CF" in calls[B]


def test_recipe_levels_of_a_multi_component_recipe(tmp_path):
    from thermochemistry_step import cbh_approach

    (tmp_path / "recipes").mkdir()
    (tmp_path / "recipes" / "toy2.yaml").write_text(
        "recipe: toy2\nversion: 1\ngeometry: {level: 'ORCA:CMP@R2SCAN-3C'}\n"
        "components:\n  a: {level: 'ORCA:QC@MP2/def2-TZVP'}\n"
        "  b: {level: 'ORCA:QC@MP2/def2-SVP'}\ncombine: |\n  E0 = a - b\n"
    )
    levels = cbh_approach.recipe_levels("SEAMM:CMP@toy2", data_path=[tmp_path])
    assert levels["energy"] == levels["recipe"] == "SEAMM:CMP@toy2"
    assert levels["components"] == {
        "a": "ORCA:QC@MP2/def2-TZVP",
        "b": "ORCA:QC@MP2/def2-SVP",
    }


def test_rerun_from_the_store_reproduces_the_composite(
    composite, monkeypatch, tmp_path
):
    """Run 1 computes every component and adds them to the store; run 2 computes
    nothing and gives the same composite record and E0 within 1e-10 E_h (the
    store keeps energies through a unit conversion, so not bitwise)."""
    import json
    import types

    import seamm_thermochemistry
    from seamm_thermochemistry import FragmentStore

    monkeypatch.setattr(
        seamm_thermochemistry, "store_is_writable", lambda path=None: (True, "test")
    )
    job = tmp_path / "job"
    job.mkdir()
    (job / "job_data.json").write_text(
        "!MolSSI job_data 1.0\n" + json.dumps({"job id": 78})
    )
    composite.flowchart = types.SimpleNamespace(root_directory=str(job))
    store = tmp_path / "fragments.db"
    P = _computed(composite, **{"fragment store": "use and add what is computed"})
    P["fragment store file"] = str(store)

    def results():
        out = {}
        for can in ("CCF", "C", "CF", "CC"):
            configuration = composite._find_fragment(composite._test_db, can)
            prefix = "#Thermochemistry#SEAMM:CMP@toy"
            out[can] = (
                _property(configuration, "composite energy" + prefix),
                _property(configuration, "composite fingerprint" + prefix),
            )
        return out

    first = composite.run_cbh(dict(P))
    computed = results()
    with FragmentStore(store) as s:
        rows = s._connect().execute("SELECT energy_fingerprint FROM energy")
        fingerprints = sorted(r[0] for r in rows.fetchall())
    # Components only, never the composite
    assert set(fingerprints) == {"fp:A/x", "fp:B/y"}

    ComponentEvaluator.calls = []
    second = composite.run_cbh(dict(P))
    assert ComponentEvaluator.calls == []
    again = results()
    for can, (E0, fingerprint) in computed.items():
        assert again[can][1] == fingerprint
        assert again[can][0] == pytest.approx(E0, abs=1e-10)
    assert second["DfH 298K"] == pytest.approx(first["DfH 298K"], abs=1e-6)
