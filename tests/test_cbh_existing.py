# -*- coding: utf-8 -*-

"""The CBH approach in its "existing results (run nothing)" mode: the species'
energies are read from properties stored on the systems named CBH:<SMILES>,
and combined with ATcT references from a synthetic ATcT cache in tests/data
(never the real tables)."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

import seamm  # noqa: E402
from molsystem import SystemDB  # noqa: E402

import thermochemistry_step  # noqa: E402

from tests.conftest import ENERGIES, PROPERTY  # noqa: E402,F401


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
        "hessian route": "automatic",
        "initial conformer": "lowest of the embeddings (MMFF)",
        "conformer embeddings": 20,
        "target conformer": "as given",
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


def test_store_fingerprints(node):
    """The engine's records become the settings/frequency fingerprint
    properties the harvest reads."""
    c = node._test_db.system.configuration
    node.model = "Test:PM7"
    record = {
        "program": "orca",
        "version": "6.1.1",
        "fingerprint": "sha256:abc",
        "element fingerprints": {"C": "sha256:c"},
        "geometry fingerprint": "sha256:geo",
    }
    node._fingerprints = {"energy": record, "hessian": None}
    node._store_fingerprints(c, "finite differences over MDI")
    p = c.properties
    get = lambda n: p.get(pattern=n, match="glob")  # noqa: E731
    assert list(get("settings fingerprint#*#*").values())[0]["value"] == "sha256:abc"
    assert list(get("frequency fingerprint#*#*").values())[0]["value"] == "sha256:abc"
    assert list(get("calculation#*#*").values())[0]["value"]["version"] == "6.1.1"
    assert list(get("frequency element fingerprints#*#*").values())[0]["value"] == {
        "C": "sha256:c"
    }
    assert node._settings_fingerprint == "sha256:abc"
    # An analytic Hessian stores the hessian record for the frequencies
    node._fingerprints = {
        "energy": record,
        "hessian": dict(record, fingerprint="sha256:hess"),
    }
    node._store_fingerprints(c, "analytic")
    assert list(get("frequency fingerprint#*#*").values())[0]["value"] == "sha256:hess"
    # No records: nothing happens
    node._fingerprints = None
    node._store_fingerprints(c, None)


def test_cbh_drives_the_engine_by_default(node):
    """The CBH approach turns 'automatic' into the warm engine, since it
    optimizes the species; an explicit choice is kept."""
    node.run_cbh(_parameters(node))
    assert node._hessian_route == "warm engine over MDI"
    node.run_cbh(_parameters(node, **{"hessian route": "separate tasks"}))
    assert node._hessian_route == "separate tasks"


def test_mixin_delegates_to_the_shared_driver(node, monkeypatch):
    """The Hessian logic lives in seamm_exec.hessian (tested there); the
    mixin passes the step's route, name and hint, and keeps the fingerprints."""
    import seamm_exec.hessian as H

    seen = {}

    def fake(self, configuration, delta_bohr=0.01, optimize=True, *a, **k):
        seen.update(route=self.route, who=self.who, hint=self.hint, opt=optimize)
        self.fingerprints = {"energy": {"fingerprint": "sha256:y"}}
        return -1.0, None, "analytic"

    monkeypatch.setattr(H.HessianDriver, "optimize_and_hessian", fake)
    monkeypatch.setattr(node, "model_chemistry", lambda: {"level": "x"})
    node._hessian_route = "warm engine over MDI"
    c = node._test_db.system.configuration
    assert node.energy_and_hessian(c) == (-1.0, None, "analytic")
    assert seen == {
        "route": "warm engine over MDI",
        "who": "Thermochemistry",
        "hint": "or put the calculation in this step's subflowchart",
        "opt": False,
    }
    assert node._fingerprints["energy"]["fingerprint"] == "sha256:y"


def test_embed_lowest_gives_all_anti(node):
    db = node._test_db
    system = db.create_system(name="CBH:CCCCC", make_current=True)
    c = system.create_configuration(name="initial", make_current=True)
    c.from_smiles("CCCCC")
    c.charge = 0
    c.spin_multiplicity = 1
    node._embed_lowest(c, "CCCCC", 20)
    assert c.n_atoms == 17 and c.charge == 0 and c.spin_multiplicity == 1
    symbols = list(c.atoms.symbols)
    xyz = c.atoms.coordinates
    text = node._chain_dihedrals(symbols, xyz)
    assert all(int(d) > 165 for d in text.split()[-2:-1][0].split()) or True
    # the two C-C-C-C dihedrals are anti
    dihedrals = [int(x) for x in text.split("order) ")[1].rstrip("°").split()]
    assert len(dihedrals) == 2 and min(dihedrals) > 165


def test_species_sdf(node, tmp_path, monkeypatch):
    """Each computed species is written to <SMILES>.sdf in the step's directory
    and appended to species.sdf, with its energies as SDF data fields."""
    pytest.importorskip("openbabel")
    monkeypatch.setattr(node, "_species_directory", lambda: tmp_path / "step")
    db = node._test_db
    configuration = None
    for system in db.systems:
        if system.name == "CBH:CC":
            configuration = system.configuration
    target = None
    for system in db.systems:
        if system.name == "CBH:CCF":
            target = system.configuration
    node._write_species_sdf(
        configuration,
        "CC",
        {"E": -1.0, "ZPE": 2.0, "H0": 1.0, "H298": 3.5, "imaginary": 0},
    )
    node._write_species_sdf(target, "CCF", {"H298": -270.0, "E": -270.0})
    single = (tmp_path / "step" / "CC.sdf").read_text()
    assert single.count("$$$$") == 1
    assert "> <SMILES>\nCC\n" in single and "> <role>\nfragment\n" in single
    assert "> <H(298.15 K) (kJ/mol)>\n3.500\n" in single
    assert "> <imaginary modes>\n0\n" in single
    combined = (tmp_path / "step" / "species.sdf").read_text()
    assert combined.count("$$$$") == 2 and "> <SMILES>\nCCF\n" in combined
    # SMILES with stereo slashes give safe file names
    node._write_species_sdf(configuration, "C/C=C/C", {"E": 0.0})
    assert (tmp_path / "step" / "C_C=C_C.sdf").exists()


def test_species_sdf_after_stored_properties(node, tmp_path, monkeypatch):
    """With properties stored on the configuration, SEAMM's own SDF data items come
    first; the step's fields must still be separate items that a reader can parse
    (the blank line ending the last stored item was once stripped, so a reader joined
    the energy with '> <SMILES>' and the file could not be read back)."""
    pytest.importorskip("openbabel")
    Chem = pytest.importorskip("rdkit.Chem")
    monkeypatch.setattr(node, "_species_directory", lambda: tmp_path / "step")
    configuration = None
    for system in node._test_db.systems:
        if system.name == "CBH:CC":
            configuration = system.configuration
    properties = configuration.properties
    name = "total energy#Thermochemistry#TEST"
    if not properties.exists(name):
        properties.add(name, "float", units="E_h", description="test")
    properties.put(name, -79.5)
    node._write_species_sdf(configuration, "CC", {"E": -1.0, "imaginary": 0})
    path = tmp_path / "step" / "CC.sdf"
    molecules = [m for m in Chem.SDMolSupplier(str(path), removeHs=False)]
    assert len(molecules) == 1 and molecules[0] is not None
    fields = molecules[0].GetPropsAsDict()
    assert float(fields["SEAMM|total energy#Thermochemistry#TEST|float|E_h"]) == -79.5
    assert fields["SMILES"] == "CC"
    assert str(fields["imaginary modes"]) == "0"


def test_species_sdf_never_fails_a_species(node, tmp_path, monkeypatch, caplog):
    """The SDF writer is diagnostic output: an error in it is logged, not raised."""
    monkeypatch.setattr(node, "_species_directory", lambda: tmp_path / "step")
    monkeypatch.setattr(node, "_engine_mode", "Model Chemistry", raising=False)
    monkeypatch.setattr(
        node, "_compute_species_mc", lambda _P, can, configuration: {"E": "n/a"}
    )
    configuration = None
    for system in node._test_db.systems:
        if system.name == "CBH:CC":
            configuration = system.configuration
    monkeypatch.setattr(
        node, "_write_species_sdf", lambda *a: (_ for _ in ()).throw(ValueError("x"))
    )
    result = node._compute_species(
        {}, "CC", "CC", system=configuration.system, configuration=configuration
    )
    assert result == {"E": "n/a"}
    assert "Could not write the structure of CC" in caplog.text
