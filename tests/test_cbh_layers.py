# -*- coding: utf-8 -*-

"""The CBH approach over the two-layer fragment store (P5c, see
seamm_thermochemistry's docs/developer_guide/campaigns/2026-11-10, §4): species
are found in the official store when the local one lacks them, an energy on an
official geometry may come from the local store (R5), and the species table
says which store each value came from."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from tests.conftest import ENERGIES  # noqa: E402
from tests.test_cbh_energy_level import (  # noqa: E402,F401
    LOW,
    OFFSET,
    FakeEvaluator,
    _computed,
    split,
)


def _put_cf(store, db, node, *, energy=True):
    from seamm_thermochemistry import EnergyRecord, GeometryRecord

    configuration = node._find_fragment(db, "CF")
    symbols = configuration.atoms.symbols
    xyz = configuration.atoms.get_coordinates(fractionals=False)
    E = ENERGIES["CF"] + OFFSET["CF"]
    g = GeometryRecord.make(
        "CF", symbols, xyz, E=E, ZPE=10.0, H=E + 15.0, model_chemistry=LOW
    )
    store.put_geometry(g)
    if energy:
        store.put_energy(
            EnergyRecord.make(
                g["geometry_hash"], ENERGIES["CF"], energy_fingerprint="fp:high"
            )
        )


@pytest.fixture()
def layers(split, monkeypatch, tmp_path):  # noqa: F811
    """A local store and an official one (version 2026.12.1) in tmp_path."""
    import seamm_thermochemistry.layered_store as L

    official = tmp_path / "fragments-official.db"
    monkeypatch.setattr(L, "official_store_path", lambda path=None: official)
    split._local = tmp_path / "fragments.db"
    split._official = official
    P = _computed(split, **{"fragment store": "use only"})
    P["fragment store file"] = str(split._local)
    split._P = P
    return split


def _official(path):
    from seamm_thermochemistry import FragmentStore

    store = FragmentStore(path, create=True, layer="official")
    store.set_meta(version="2026.12.1")
    return store


def test_a_species_from_the_official_store(layers):
    with _official(layers._official) as store:
        _put_cf(store, layers._test_db, layers)
    FakeEvaluator.calls = []
    data = layers.run_cbh(dict(layers._P))
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert "CF" not in layers._computed
    assert all("CBH:CF" not in call for call in FakeEvaluator.calls)
    assert layers._store_layer["CF"] == "official"
    row = {"layer": "official", "verified": False}
    assert layers._store_source(layers._fragment_store(layers._P), row) == (
        "store (official 2026.12.1, unverified)"
    )


def test_a_local_energy_on_an_official_geometry(layers):
    """R5: the geometry from the official store, its energy-level energy from the
    local one (same geometry hash): used, and both stores named."""
    from seamm_thermochemistry import FragmentStore

    with _official(layers._official) as store:
        _put_cf(store, layers._test_db, layers, energy=False)
    with FragmentStore(layers._local, create=True) as store:
        _put_cf(store, layers._test_db, layers)  # the user's energy on it
    # Only the local energy row is wanted: drop the local geometry row
    with FragmentStore(layers._local, create=True) as store:
        store._connect(create=True).execute("DELETE FROM geometry")
        store._connect(create=True).commit()
    seen = {}
    real = layers._energy_level_energies

    def spy(mc, store, target, order, configurations, energies, sources):
        real(mc, store, target, order, configurations, energies, sources)
        seen.update(sources)

    layers._energy_level_energies = spy
    FakeEvaluator.calls = []
    data = layers.run_cbh(dict(layers._P))
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert all("CBH:CF" not in call for call in FakeEvaluator.calls)
    assert "geometry official 2026.12.1, energy local" in seen["CF"]


def test_local_is_preferred_and_without_an_official_store_nothing_changes(layers):
    from seamm_thermochemistry import FragmentStore

    with FragmentStore(layers._local, create=True) as store:
        _put_cf(store, layers._test_db, layers)
    FakeEvaluator.calls = []
    layers.run_cbh(dict(layers._P))
    assert layers._store_layer["CF"] == "local"
    store = layers._fragment_store(layers._P)
    assert store.official is None
    assert layers._store_source(store, {"layer": "local", "verified": True}) == (
        "store"
    )
