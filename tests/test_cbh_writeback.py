# -*- coding: utf-8 -*-

"""Write-back (P5b of the fragment store, see seamm_thermochemistry's
docs/developer_guide/campaigns/2026-11-10): with "use and add what is computed"
the species the step computes go into the fragment store at once, so later
iterations and jobs find them; when the store cannot be written here, the job is
marked to be harvested when it returns."""

import json
import types

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from tests.test_cbh_energy_level import (  # noqa: E402,F401
    FakeEvaluator,
    _computed,
    split,
)


@pytest.fixture()
def adding(split, monkeypatch, tmp_path):  # noqa: F811
    """The split node adding to a store in tmp_path, which is writable, in a job
    whose directory is tmp_path/job."""
    import seamm_thermochemistry

    # The step imports the name from the package
    monkeypatch.setattr(
        seamm_thermochemistry, "store_is_writable", lambda path=None: (True, "test")
    )
    job = tmp_path / "job"
    job.mkdir()
    (job / "job_data.json").write_text(
        "!MolSSI job_data 1.0\n" + json.dumps({"job id": 77})
    )
    split.flowchart = types.SimpleNamespace(root_directory=str(job))
    store = tmp_path / "fragments.db"
    P = _computed(split, **{"fragment store": "use and add what is computed"})
    P["fragment store file"] = str(store)
    split._P = P
    split._store = store
    split._job = job
    return split


def test_computed_species_are_added_and_reused(adding):
    from seamm_thermochemistry import FragmentStore

    FakeEvaluator.calls = []
    adding.run_cbh(dict(adding._P))
    with FragmentStore(adding._store) as store:
        rows = store._connect().execute("SELECT job, energy_fingerprint FROM energy")
        rows = rows.fetchall()
    # The energy level's single points of the target and the three fragments
    assert len(rows) == 4
    assert {r[1] for r in rows} == {"fp:high"}
    assert all(r[0].startswith("77:") for r in rows)
    # A second job finds them: only the target runs at the energy level
    FakeEvaluator.calls = []
    adding.run_cbh(dict(adding._P))
    assert FakeEvaluator.calls == [["CBH:CCF"]]


def test_not_writable_marks_the_job_for_harvest(adding, monkeypatch):
    import seamm_thermochemistry

    monkeypatch.setattr(
        seamm_thermochemistry,
        "store_is_writable",
        lambda path=None: (False, "on a gpfs filesystem"),
    )
    adding.run_cbh(dict(adding._P))
    assert not adding._store.exists()
    marker = json.loads((adding._job / "harvest_fragments.json").read_text())
    assert sorted(marker["species"]) == ["C", "CC", "CCF", "CF"]
    assert "gpfs" in marker["reason"]


@pytest.mark.parametrize("choice", ["use only", "use if available", "ignore"])
def test_other_choices_add_nothing(adding, choice):
    P = dict(adding._P, **{"fragment store": choice})
    adding.run_cbh(P)
    assert not adding._store.exists()
    assert not (adding._job / "harvest_fragments.json").exists()


def test_write_back_never_fails_the_step(adding, monkeypatch):
    import seamm_thermochemistry

    def broken(*args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(seamm_thermochemistry, "harvest_configuration", broken)
    data = adding.run_cbh(dict(adding._P))
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    marker = json.loads((adding._job / "harvest_fragments.json").read_text())
    assert "disk full" in marker["reason"]


def test_the_default_is_to_add():
    from thermochemistry_step.thermochemistry_parameters import (
        ThermochemistryParameters,
    )

    p = ThermochemistryParameters.parameters["fragment store"]
    assert p["default"] == "use and add what is computed"
    assert p["enumeration"] == ("use and add what is computed", "use only", "ignore")
