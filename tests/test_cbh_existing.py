# -*- coding: utf-8 -*-

"""The CBH approach in its "existing results (run nothing)" mode: the species'
energies are read from properties stored on the systems named CBH:<SMILES>,
and combined with ATcT references from a synthetic ATcT cache in tests/data
(never the real tables)."""

from pathlib import Path

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

import seamm  # noqa: E402
from molsystem import SystemDB  # noqa: E402
from seamm_thermochemistry import atct  # noqa: E402

import thermochemistry_step  # noqa: E402

FIXTURE = Path(__file__).parent / "data" / "atct_fixture.json"
PROPERTY = "total energy#Test#PM7"

# Synthetic energies, kJ/mol, taken as enthalpies at 298.15 K (no thermal
# corrections). Fluoroethane CBH-1: CCF + C -> CF + CC, so
# DrH = E(CF) + E(CC) - E(CCF) - E(C) = -230 - 80 + 262 + 50 = 2.0
ENERGIES = {"CCF": -262.0, "C": -50.0, "CF": -230.0, "CC": -80.0}


@pytest.fixture()
def atct_base(tmp_path):
    base = tmp_path / "ATcT"
    atct.download_all(base, source=FIXTURE)
    return base


@pytest.fixture()
def node(tmp_path, atct_base, monkeypatch):
    db = SystemDB(filename=f"file:{tmp_path / 'seamm.db'}")
    target = None
    for smiles, energy in ENERGIES.items():
        system = db.create_system(name=f"CBH:{smiles}", make_current=True)
        configuration = system.create_configuration(name="initial", make_current=True)
        configuration.from_smiles(smiles)
        configuration.properties.add(PROPERTY, "float", units="kJ/mol", noerror=True)
        configuration.properties.put(PROPERTY, energy)
        if smiles == "CCF":
            target = (system, configuration)
    db.system = target[0]

    node = thermochemistry_step.Thermochemistry()
    variables = {"_system_db": db}
    monkeypatch.setattr(node, "get_variable", lambda name: variables[name])
    monkeypatch.setattr(node, "get_system_configuration", lambda *a, **k: target)
    monkeypatch.setattr(node, "get_table", lambda name: seamm.Table.create(db, name))
    monkeypatch.setattr(node, "_cite_cbh", lambda tables: None)
    node._test_db = db
    node._test_base = atct_base
    yield node
    db.close()


def _parameters(node, **changes):
    P = {
        "cbh rungs": "1",
        "protect": "default",
        "cbh energies": "existing results (run nothing)",
        "thermal corrections": "none (the energy is the enthalpy, e.g. MOPAC)",
        "energy property": "total energy#*#*",
        "ZPE property": "ZPE#*#*",
        "enthalpy property": "H(298.15 K)#*#*",
        "low-level energy property": "",
        "fragment store": "ignore",
        "fragment store file": "default",
        "reference values": "ATcT",
        "ladder": "",
        "ATcT uncertainty threshold": 0.5,
        "ATcT version": "9.999",
        "ATcT directory": str(node._test_base),
        "covariances": "cached only",
        "compare with ATcT": "yes",
    }
    P.update(changes)
    return P


def test_existing_results(node):
    data = node.run_cbh(_parameters(node))
    assert data["CBH rung"] == 1
    assert data["CBH reaction"] == "CBH-1: C + CCF -> CC + CF"
    # ATcT 298 K references of the fixture: CF -235.55 + CC -83.91 - C -74.513
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert "DfH 0K" not in data  # no 0 K value without thermal corrections
    assert data["DfH uncertainty"] == pytest.approx(
        (0.29**2 + 0.14**2 + 0.043**2) ** 0.5, abs=1e-6
    )
    assert data["ATcT DfH 298K"] == pytest.approx(-272.11)
    assert data["ATcT version"] == "9.999"
    table = seamm.Table(node._test_db, "CBH")
    assert table.n_rows == 1
    assert data["CBH results"]["1"]["uncertainty method"].startswith("quadrature")


# The n-propyl radical, CBH-1: [CH2]CC + C -> [CH2]C + CC (the ethyl
# radical's own CBH-1 is trivial), DrH = E([CH2]C) + E(CC) - E([CH2]CC) - E(C)
# = 130 - 80 - 110 + 50 = -10
RADICAL_ENERGIES = {"[CH2]CC": 110.0, "C": -50.0, "[CH2]C": 130.0, "CC": -80.0}


@pytest.fixture()
def radical_node(tmp_path, atct_base, monkeypatch):
    db = SystemDB(filename=f"file:{tmp_path / 'radical.db'}")
    target = None
    for smiles, energy in RADICAL_ENERGIES.items():
        system = db.create_system(name=f"CBH:{smiles}", make_current=True)
        configuration = system.create_configuration(name="initial", make_current=True)
        configuration.from_smiles(smiles)
        configuration.properties.add(PROPERTY, "float", units="kJ/mol", noerror=True)
        configuration.properties.put(PROPERTY, energy)
        if smiles == "[CH2]CC":
            target = (system, configuration)
    db.system = target[0]
    # from_smiles sets the multiplicity of a radical correctly ...
    assert target[1].spin_multiplicity == 2 and target[1].charge == 0

    node = thermochemistry_step.Thermochemistry()
    variables = {"_system_db": db}
    monkeypatch.setattr(node, "get_variable", lambda name: variables[name])
    monkeypatch.setattr(node, "get_system_configuration", lambda *a, **k: target)
    monkeypatch.setattr(node, "get_table", lambda name: seamm.Table.create(db, name))
    monkeypatch.setattr(node, "_cite_cbh", lambda tables: None)
    node._test_db = db
    node._test_base = atct_base
    node._test_target = target[1]
    yield node
    db.close()


def test_radical_target(radical_node, monkeypatch):
    node = radical_node
    # ... but to_smiles of a radical gives an anion (molsystem, 2026-10-28),
    # so stand in for the fixed conversion here
    monkeypatch.setattr(node._test_target, "to_smiles", lambda **k: "[CH2]CC")
    data = node.run_cbh(_parameters(node))
    from seamm_thermochemistry import cbh

    reaction = cbh.cbh_reactions("[CH2]CC", rungs=[1])[1]
    assert dict(reaction.references) == {"[CH2]C": 1, "CC": 1, "C": -1}
    assert data["CBH reaction"] == str(reaction)
    # 298 K references: [CH2]C 120.75 + CC -83.91 - C -74.513 = 111.353
    assert data["DfH 298K"] == pytest.approx(111.353 + 10.0, abs=1e-6)
    assert "ATcT DfH 298K" not in data  # the propyl radical is not in the fixture


def test_radical_smiles_guard(radical_node, monkeypatch):
    """A SMILES that disagrees with the configuration's charge and spin (as
    molsystem before 2026.10.28 gave for a radical, [CH2-]CC) must stop the
    step rather than build the reactions of the anion."""
    monkeypatch.setattr(radical_node._test_target, "to_smiles", lambda **k: "[CH2-]CC")
    with pytest.raises(RuntimeError, match=r"implies charge -1 and 0 unpaired"):
        radical_node.run_cbh(_parameters(radical_node))


def test_radical_target_from_structure(radical_node):
    """With molsystem's radical-aware SMILES (2026.10.28), the target's SMILES
    comes straight from the structure."""
    smiles = radical_node._test_target.to_smiles(canonical=True, isomeric=True)
    if "-" in smiles:
        pytest.skip(f"molsystem gives {smiles} for the n-propyl radical")
    data = radical_node.run_cbh(_parameters(radical_node))
    assert data["DfH 298K"] == pytest.approx(111.353 + 10.0, abs=1e-6)


def test_fragment_charge_and_multiplicity():
    from thermochemistry_step.cbh_approach import cbh_charge, cbh_multiplicity

    assert (cbh_charge("[CH3]"), cbh_multiplicity("[CH3]")) == (0, 2)
    assert (cbh_charge("[NH3+]C"), cbh_multiplicity("[NH3+]C")) == (1, 1)
    assert (cbh_charge("CC(=O)[O-]"), cbh_multiplicity("CC(=O)[O-]")) == (-1, 1)
    assert (cbh_charge("C[N+](=O)[O-]"), cbh_multiplicity("C[N+](=O)[O-]")) == (0, 1)
    assert (cbh_charge("[CH2]"), cbh_multiplicity("[CH2]")) == (0, 3)


def test_missing_fragment(node):
    with pytest.raises(RuntimeError, match="No property matching"):
        node.run_cbh(_parameters(node, **{"energy property": "nothing#*#*"}))


def test_fragment_from_store(node, monkeypatch, tmp_path):
    """A fragment with no system is served from the fragment store at the
    step's model chemistry; an unverified record is reported as such."""
    from seamm_thermochemistry import FragmentRecord, FragmentStore

    path = tmp_path / "fragments.db"
    with FragmentStore(path, create=True) as store:
        store.put(FragmentRecord.make("CF", E=-230.0, model_chemistry="Test:PM7"))
    find = node._find_fragment
    monkeypatch.setattr(
        node,
        "_find_fragment",
        lambda db, can: None if can == "CF" else find(db, can),
    )
    node.model = "Test:PM7"
    P = _parameters(node, **{"fragment store": "use if available"})
    P["fragment store file"] = str(path)
    data = node.run_cbh(P)
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    # At another model chemistry the store has nothing
    node.model = "Other:PM6"
    with pytest.raises(RuntimeError, match="nothing in the fragment store"):
        node.run_cbh(P)


def _ladder(tmp_path, atct_base, shift):
    """A ladder whose CH3F and C2H6 values are ATcT's shifted by `shift`."""
    from seamm_thermochemistry import ATcT, build_ladder, FragmentRecord, FragmentStore

    tables = ATcT(base=atct_base)
    store_path = tmp_path / "fragments.db"
    level = "High:F12"
    with FragmentStore(store_path, create=True) as store:
        for smiles in ("[H][H]", "C", "F", "CC", "CF"):
            h0 = tables.lookup(smiles)["dfH0_0K"] + (
                shift if smiles in ("CC", "CF") else 0
            )
            store.put(
                FragmentRecord.make(
                    smiles, E=h0 - 10.0, ZPE=10.0, H=h0 + 5.0, model_chemistry=level
                )
            )
    ladder = build_ladder(
        ["CC", "CF"],
        name="test",
        energy_level=level,
        zpe_scale=1.0,
        atct_base=atct_base,
        store_path=store_path,
    )
    ladder.save(store_path)
    return store_path


def test_laddered_references(node, tmp_path):
    """Fluoroethane CBH-1 with laddered CH3F and C2H6 (both 1 kJ/mol above
    ATcT): DfH shifts by +2, the ATcT part of the uncertainty is that of the
    anchors (2 CH4 + HF - H2 ... ) and the method term is added."""
    store_path = _ladder(tmp_path, node._test_base, 1.0)
    P = _parameters(
        node,
        **{
            "reference values": "laddered",
            "ladder": "test",
            "fragment store file": str(store_path),
        },
    )
    P["thermal corrections"] = "harmonic"
    P["ZPE property"] = "none"
    P["enthalpy property"] = "none"
    # With no ZPE/H properties only E is known: use 'none' corrections instead
    P["thermal corrections"] = "none (the energy is the enthalpy, e.g. MOPAC)"
    data = node.run_cbh(P)
    # ATcT: CF -235.55 + CC -83.91 - C -74.513 = -244.947; laddered at 298 K:
    # the ladder shifts the 0 K values by 1; at 298 K H = h0 + 5 for every
    # species, and the shift carries through
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0 + 2.0, abs=0.05)
    r = data["CBH results"]["1"]
    assert "ladder method term" in r["uncertainty method"]
    assert r["uncertainty ladder method"] == pytest.approx(1.0 * 2**0.5, abs=1e-6)
    # ATcT part: CH3F -> CH4 + HF - H2 anchors; C2H6 -> 2 CH4 - H2; minus CH4
    # = 2 CH4 + HF - 2 H2: sqrt((2*0.043)^2 + 0.018^2) (H2 exact)
    assert r["uncertainty ATcT"] == pytest.approx(
        ((2 * 0.043) ** 2 + 0.018**2) ** 0.5, abs=1e-6
    )


def test_threshold_policy(node, tmp_path):
    store_path = _ladder(tmp_path, node._test_base, 1.0)
    P = _parameters(
        node,
        **{
            "reference values": "ATcT below the threshold, else laddered",
            "ladder": "test",
            "fragment store file": str(store_path),
            # CC (0.14) from ATcT, CF (0.29) laddered
            "ATcT uncertainty threshold": 0.2,
        },
    )
    data = node.run_cbh(P)
    r = data["CBH results"]["1"]
    assert r["uncertainty ladder method"] == pytest.approx(1.0, abs=1e-6)


def test_missing_fragment_system(node, monkeypatch):
    find = node._find_fragment
    monkeypatch.setattr(
        node,
        "_find_fragment",
        lambda db, can: None if can == "CF" else find(db, can),
    )
    with pytest.raises(RuntimeError, match="CBH:CF"):
        node.run_cbh(_parameters(node))
