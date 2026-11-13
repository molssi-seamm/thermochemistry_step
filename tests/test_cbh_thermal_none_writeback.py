# -*- coding: utf-8 -*-

"""Write-back without thermal corrections ("the energy is the enthalpy", e.g.
MOPAC's enthalpy of formation): each computed species keeps its energy as the
property the fragment store's harvest reads, so it is stored and reused, and a
species of which nothing could be stored is reported with the reason."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from seamm_util import Q_  # noqa: E402

MODEL = "MOPAC:SQM@PM7"
NONE = "none (the energy is the enthalpy, e.g. MOPAC)"
#: Methane, C-H 1.089 Å
CH4 = [
    ("C", (0.0, 0.0, 0.0)),
    ("H", (0.629, 0.629, 0.629)),
    ("H", (-0.629, -0.629, 0.629)),
    ("H", (-0.629, 0.629, -0.629)),
    ("H", (0.629, -0.629, -0.629)),
]


@pytest.fixture()
def methane():
    from molsystem import SystemDB

    db = SystemDB(filename="file:tnone?mode=memory&cache=shared")
    configuration = db.create_system().create_configuration()
    configuration.atoms.append(
        x=[p[0] for _, p in CH4],
        y=[p[1] for _, p in CH4],
        z=[p[2] for _, p in CH4],
        symbol=[s for s, _ in CH4],
    )
    yield configuration
    db.close()


def _species_mc(node, configuration, monkeypatch, E=-49.6, fingerprint="mopac:pm7"):
    """The Model Chemistry route for one species without thermal corrections,
    with the engine stubbed: optimized, energy E (kJ/mol), fingerprint."""
    monkeypatch.setattr(node, "optimize_and_energy", lambda *a, **k: E)
    node.model = MODEL
    node._fingerprints = (
        {"energy": {"fingerprint": fingerprint, "program": "MOPAC"}}
        if fingerprint
        else None
    )
    P = {
        "thermal corrections": NONE,
        "optimize fragments": "yes",
        "fragment convergence": Q_(0.05, "eV/Å"),
        "fragment max steps": 100,
    }
    return P, node._compute_species_mc(P, "C", configuration)


def test_the_energy_is_stored_for_the_harvest(node, methane, monkeypatch):
    P, result = _species_mc(node, methane, monkeypatch)
    assert result == {"H298": -49.6, "E": -49.6}
    name = f"total energy#Thermochemistry#{MODEL}"
    value = methane.properties.get(name)[name]["value"]
    assert Q_(value, "E_h").m_as("kJ/mol") == pytest.approx(-49.6)
    # Nothing thermal: no harmonic run may take it for one
    assert not methane.properties.get(pattern="ZPE#*", match="glob")
    assert not methane.properties.get(pattern="H(298.15 K)#*", match="glob")


def test_harvested_and_found_again(node, methane, monkeypatch, tmp_path):
    from seamm_thermochemistry import FragmentStore, harvest_configuration
    from seamm_thermochemistry.fragment_store import new_counts

    P, _ = _species_mc(node, methane, monkeypatch)
    counts = new_counts()
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        harvest_configuration(
            store,
            methane,
            "C",
            charge=0,
            multiplicity=1,
            energy="total energy#Thermochemistry#*",
            zpe="ZPE#Thermochemistry#*",
            enthalpy="H(298.15 K)#Thermochemistry#*",
            fingerprint="settings fingerprint#Thermochemistry#*",
            frequency_fingerprint="frequency fingerprint#Thermochemistry#*",
            write_back=True,
            counts=counts,
        )
        assert counts["inserted"] == 1 and counts["energy"] == {"inserted": 1}
        assert counts["geometry"] == {}, "no geometry row without thermal data"
        # The next species or job finds it, as the enthalpy
        monkeypatch.setattr(
            node, "_geometry_level_fingerprints", lambda can: (["mopac:pm7"], None)
        )
        found, source = node._energies_from_store(P, store, "C")
        assert found == {"E": pytest.approx(-49.6), "H298": pytest.approx(-49.6)}
        # A harmonic run of the same level does not take it: no H
        harmonic = dict(P, **{"thermal corrections": "harmonic"})
        assert node._energies_from_store(harmonic, store, "C") is None


def test_unfingerprinted_engines_are_stored_unverified(
    node, methane, monkeypatch, tmp_path
):
    from seamm_thermochemistry import FragmentStore, harvest_configuration

    P, _ = _species_mc(node, methane, monkeypatch, fingerprint=None)
    with FragmentStore(tmp_path / "f.db", create=True) as store:
        counts = harvest_configuration(
            store, methane, "C", energy="total energy#Thermochemistry#*"
        )
        assert counts["inserted"] == 1
        record = store.get("C", model_chemistry=MODEL)
        assert not record["verified"] and record["E"] == pytest.approx(-49.6)


def test_a_species_not_added_is_reported(node, methane, monkeypatch, tmp_path):
    """Write-back says why a computed species was not stored, instead of the
    bare 'energies: none'."""
    import seamm_thermochemistry

    monkeypatch.setattr(
        seamm_thermochemistry, "store_is_writable", lambda path=None: (True, "test")
    )
    monkeypatch.setattr(node, "_job_root", lambda: None)
    printed = []
    monkeypatch.setattr(
        "thermochemistry_step.cbh_approach.printer.important",
        lambda text: printed.append(str(text)),
    )
    P = {
        "fragment store": "use and add what is computed",
        "fragment store file": str(tmp_path / "f.db"),
    }
    node.model = MODEL
    node._write_back(P, {"C": methane})  # no energy property on it
    text = "\n".join(printed)
    assert "species: none" in text and "energies: none" in text
    assert "C was not added: no property matching 'total energy#" in text

    # With the energy: added, and nothing reported missing
    printed.clear()
    _species_mc(node, methane, monkeypatch)
    node._write_back(P, {"C": methane})
    text = "\n".join(printed)
    assert "species: 1 inserted" in text and "energies: 1 inserted" in text
    assert "was not added" not in text
