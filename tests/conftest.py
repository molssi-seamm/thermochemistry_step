# -*- coding: utf-8 -*-

"""Shared fixtures: a synthetic ATcT cache and a Thermochemistry node over a
system database with the fluoroethane CBH-1 species."""

from pathlib import Path

import pytest

FIXTURE = Path(__file__).parent / "data" / "atct_fixture.json"
PROPERTY = "total energy#Test#PM7"

# Synthetic energies, kJ/mol, taken as enthalpies at 298.15 K (no thermal
# corrections). Fluoroethane CBH-1: CCF + C -> CF + CC, so
# DrH = E(CF) + E(CC) - E(CCF) - E(C) = -230 - 80 + 262 + 50 = 2.0
ENERGIES = {"CCF": -262.0, "C": -50.0, "CF": -230.0, "CC": -80.0}


@pytest.fixture()
def atct_base(tmp_path):
    from seamm_thermochemistry import atct

    base = tmp_path / "ATcT"
    atct.download_all(base, source=FIXTURE)
    return base


@pytest.fixture()
def node(tmp_path, atct_base, monkeypatch):
    pytest.importorskip("rdkit")
    import seamm
    from molsystem import SystemDB

    import thermochemistry_step

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
