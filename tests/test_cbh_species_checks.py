# -*- coding: utf-8 -*-

"""Checks on each computed species (PPE1694 jobs 588/594, 2026-11-13): an
optimization that does not converge is continued once, then flagged and kept out
of the fragment store; a target whose bonds change while it is optimized fails
(dimethyl sulfate became SO2 + CH3OOCH3 at PM6-ORG); a fragment that does so is
flagged and kept out of the store."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from tests.test_cbh_writeback import adding  # noqa: E402,F401
from tests.test_cbh_energy_level import (  # noqa: E402,F401
    ENERGIES,
    FakeEvaluator,
    OFFSET,
    _computed,
    split,
)


def test_an_unconverged_optimization_is_continued(node):
    calls = []

    def run(converged):
        def once():
            calls.append(1)
            node._last_optimization = {"converged": converged.pop(0)}
            return len(calls)

        return once

    # Converges when continued: the second run's result, converged
    assert node._run_optimization("CCO", run([False, True]), True, 200) == (2, True)
    # Still not converged after continuing: flagged
    calls.clear()
    assert node._run_optimization("CCO", run([False, False]), True, 200) == (2, False)
    # Converged at once: one run
    calls.clear()
    assert node._run_optimization("CCO", run([True]), True, 200) == (1, True)
    # Nothing optimized, or a driver that does not say: unknown, one run
    calls.clear()
    assert node._run_optimization("CCO", run([True]), False, 200) == (1, None)

    def silent():
        node._last_optimization = None
        return 7

    assert node._run_optimization("CCO", silent, True, 200) == (7, None)


@pytest.fixture()
def ethanol():
    from molsystem import SystemDB

    db = SystemDB(filename="file:checks?mode=memory&cache=shared")

    def make(name, stretch=0.0):
        system = db.create_system(name=name)
        configuration = system.create_configuration()
        configuration.from_smiles("CCO")
        if stretch:
            xyz = configuration.atoms.get_coordinates(fractionals=False)
            symbols = configuration.atoms.symbols
            o = symbols.index("O")
            xyz[o] = [xyz[o][0] + stretch, xyz[o][1], xyz[o][2]]
            configuration.atoms.set_coordinates(xyz, fractionals=False)
        return configuration

    yield make
    db.close()


def test_a_target_that_changed_fails(node, ethanol):
    result = {"E": 1.0}
    node._check_species("CCO", ethanol("ethanol"), result)  # intact: nothing
    assert "broken" not in result
    with pytest.raises(RuntimeError, match="changed while it was optimized"):
        node._check_species("CCO", ethanol("ethanol", stretch=2.5), {"E": 1.0})


def test_a_fragment_that_changed_is_flagged(node, ethanol):
    result = {"E": 1.0}
    node._check_species("CCO", ethanol("CBH:CCO", stretch=2.5), result)
    assert result["broken"]
    assert node._storable("CCO", result) is False
    assert node._storable("CCO", {"E": 1.0, "converged": False}) is False
    assert node._storable("CCO", {"E": 1.0, "converged": True}) is True
    flags = node._flags(
        {
            "CCO": result,
            "C": {"converged": False},
            "CC": {"converged": True},
            "CF": {},
        }
    )
    assert flags == {
        "CBH unconverged species": "C",
        "CBH changed fragments": "CCO",
        "CBH saddle-point species": "",
    }


def test_unconverged_species_are_flagged_and_not_stored(
    adding, monkeypatch  # noqa: F811
):
    """A fragment whose optimization did not converge is used, flagged in the
    results, and not written back."""
    from seamm_thermochemistry import FragmentStore

    def compute(P, can, generated, system=None, configuration=None, embed=False):
        E = ENERGIES[can] + OFFSET[can]
        result = {"E": E, "ZPE": 10.0, "H0": E + 10.0, "H298": E + 15.0}
        if can == "CF":
            result["converged"] = False
        return result

    monkeypatch.setattr(adding, "_compute_species", compute)
    FakeEvaluator.calls = []
    data = adding.run_cbh(dict(adding._P))
    assert data["CBH unconverged species"] == "CF"
    assert data["CBH changed fragments"] == ""
    with FragmentStore(adding._store) as store:
        species = {
            r[0] for r in store._connect().execute("SELECT species FROM fragment")
        }
    assert "CF" not in species and {"C", "CC", "CCF"} <= species
