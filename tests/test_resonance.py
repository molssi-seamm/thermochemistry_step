# -*- coding: utf-8 -*-

"""Matching a target to ATcT across resonance forms (SiO, OBrO), and the
structure's spin multiplicity over the SMILES' radical count (SiCl2, BH)."""

import json

import pytest

pytest.importorskip("rdkit")

from thermochemistry_step.resonance import (  # noqa: E402
    lookup_target,
    resonance_key,
)


def _entry(atct_id, name, formula, smiles, dfh0, dfh298, unc="0.10", charge="0"):
    return {
        "ATcT_TN_Version": "9.998",
        "ATcT_ID": atct_id,
        "Name": name,
        "Formula": f"{formula}  (g)",
        "Delta_Hf_0K": dfh0,
        "Delta_Hf_298K": dfh298,
        "Delta_Hf298K_uncertainty": unc,
        "unit": "kJ/mol",
        "mass": "1.0",
        "mass_uncertainty": "0.0",
        "SMILES": smiles,
        "CASRN": atct_id.split("*")[0],
        "charge": charge,
    }


# The SMILES as ATcT 1.228 writes them; the values are made up
ENTRIES = [
    _entry("10097-28-6*0", "Silicon monoxide", "SiO", "[Si]=O", "-100.0", "-99.0"),
    _entry("21255-83-4*0", "Bromine dioxide", "BrO2", "[O]Br[O]", "160.0", "155.0"),
    _entry("630-08-0*0", "Carbon monoxide", "CO", "[C-]#[O+]", "-113.8", "-110.5"),
    _entry("74-90-8*0", "Hydrogen cyanide", "HCN", "C#N", "129.7", "129.6"),
    _entry("6914-07-4*0", "Hydrogen isocyanide", "HNC", "[C]=N", "192.0", "191.5"),
    _entry("74-82-8*0", "Methane", "CH4", "C", "-66.550", "-74.513", "0.043"),
]


@pytest.fixture()
def tables(tmp_path):
    from seamm_thermochemistry import atct

    source = tmp_path / "atct.json"
    source.write_text(json.dumps(ENTRIES))
    base = tmp_path / "ATcT"
    atct.download_all(base, source=source)
    return atct.ATcT(base=base)


@pytest.mark.parametrize(
    "a, b",
    [
        ("[O+]#[Si-]", "[Si]=O"),
        ("[O][Br+][O-]", "[O]Br[O]"),
        ("[C-]#[O+]", "[C]=O"),
        ("[C-]#[NH+]", "[C]=N"),
    ],
)
def test_resonance_forms_share_a_key(a, b):
    assert resonance_key(a) == resonance_key(b)


@pytest.mark.parametrize(
    "a, b",
    [
        ("C#N", "[C]=N"),  # HCN and HNC are isomers
        ("CO", "C=O"),  # different hydrogen counts
        ("[NH4+]", "N"),  # different charges
    ],
)
def test_different_species_differ(a, b):
    assert resonance_key(a) != resonance_key(b)


@pytest.mark.parametrize(
    "target, multiplicity, atct_id, via",
    [
        ("[O+]#[Si-]", 1, "10097-28-6*0", "O=[Si]"),
        ("[O][Br+][O-]", 2, "21255-83-4*0", "[O]Br[O]"),
        ("[C-]#[O+]", 1, "630-08-0*0", None),  # a direct match
        ("C#N", 1, "74-90-8*0", None),  # HCN, not HNC
        ("[C-]#[NH+]", 1, "6914-07-4*0", "[C]=N"),  # HNC, not HCN
    ],
)
def test_lookup_target(tables, target, multiplicity, atct_id, via):
    entry, matched = lookup_target(tables, target, 0, multiplicity)
    assert entry["id"] == atct_id
    assert matched == via


def test_no_match(tables):
    assert lookup_target(tables, "CC", 0, 1) == (None, None)
    # The charge is part of the species
    assert lookup_target(tables, "[O+]#[Si-]", 1, 2) == (None, None)


def test_several_matches_are_not_used(tables, tmp_path):
    from seamm_thermochemistry import atct

    # Two ATcT species with one resonance-invariant key: no match is chosen
    entries = ENTRIES + [
        _entry("99999-99-9*0", "Silicon monoxide, ionic", "SiO", "[Si+][O-]", "0", "0")
    ]
    source = tmp_path / "atct2.json"
    source.write_text(json.dumps(entries))
    base = tmp_path / "ATcT2"
    atct.download_all(base, source=source)
    tables2 = atct.ATcT(base=base)
    assert lookup_target(tables2, "[O+]#[Si-]", 0, 1) == (None, None)


# -- the structure's multiplicity over the SMILES' radicals --------------------


class _Configuration:
    def __init__(self, smiles):
        self.smiles = smiles

    def to_smiles(self, **kwargs):
        return self.smiles


@pytest.mark.parametrize("smiles", ["Cl[Si]Cl", "[BH]"])
def test_singlet_written_as_diradical(smiles):
    from thermochemistry_step.cbh_approach import CBHApproach

    assert CBHApproach._target_smiles(_Configuration(smiles), 0, 1) == smiles


def test_parity_and_charge_still_checked():
    from thermochemistry_step.cbh_approach import CBHApproach

    # A doublet cannot be a singlet: the electron count disagrees
    with pytest.raises(RuntimeError, match="implies charge 0 and 1 unpaired"):
        CBHApproach._target_smiles(_Configuration("[CH3]"), 0, 1)
    # A radical rendered as an anion
    with pytest.raises(RuntimeError, match="implies charge -1"):
        CBHApproach._target_smiles(_Configuration("[CH2-]C"), 0, 2)


# -- end to end: a target that ATcT writes as another resonance form -----------


def test_sio_target_is_its_atct_value(tables, tmp_path, monkeypatch):
    """SiO from a structure is [O+]#[Si-]; ATcT writes [Si]=O. Every rung is
    trivial, so the result is the ATcT value, found through the resonance form."""
    import seamm
    from molsystem import SystemDB

    import thermochemistry_step
    from tests.test_cbh_existing import _parameters

    db = SystemDB(filename=f"file:{tmp_path / 'sio.db'}")
    system = db.create_system(name="SiO", make_current=True)
    configuration = system.create_configuration(name="initial", make_current=True)
    configuration.from_smiles("[O+]#[Si-]")
    monkeypatch.setattr(configuration, "to_smiles", lambda **k: "[O+]#[Si-]")
    db.system = system

    node = thermochemistry_step.Thermochemistry()
    variables = {"_system_db": db}
    monkeypatch.setattr(node, "get_variable", lambda name: variables[name])
    monkeypatch.setattr(
        node, "get_system_configuration", lambda *a, **k: (system, configuration)
    )
    monkeypatch.setattr(node, "get_table", lambda name: seamm.Table.create(db, name))
    monkeypatch.setattr(node, "_cite_cbh", lambda tables: None)
    node._test_base = tables.directory.parent

    try:
        data = node.run_cbh(_parameters(node, **{"ATcT version": "9.998"}))
    finally:
        db.close()
    assert data["CBH reaction"] == "target is a reference species (ATcT 10097-28-6*0)"
    assert data["DfH 298K"] == pytest.approx(-99.0)
    assert data["DfH 0K"] == pytest.approx(-100.0)
