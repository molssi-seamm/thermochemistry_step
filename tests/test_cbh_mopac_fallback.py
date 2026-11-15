# -*- coding: utf-8 -*-

"""Decisions of 2026-11-15 (Paul): a target whose optimization did not converge
fails, like one whose bonds changed; and with MOPAC, when no CBH rung of 2 or
higher is usable, the result is MOPAC's own enthalpy of formation, reported as
"MOPAC ΔHf (no CBH rung ≥2)" (below CBH-2, CBH with MOPAC is worse than MOPAC's
own value: PPE1694, PM7 CBH-0 353, CBH-1 35, own 13-24 kJ/mol)."""

from types import SimpleNamespace

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from thermochemistry_step import cbh_approach  # noqa: E402

from tests.conftest import ENERGIES  # noqa: E402
from tests.test_cbh_existing import _parameters  # noqa: E402
from tests.test_cbh_species_checks import ethanol  # noqa: E402,F401

MOPAC = "MOPAC:SQM@PM7"


def test_an_unconverged_target_fails(node, ethanol):  # noqa: F811
    with pytest.raises(RuntimeError, match="did not converge"):
        node._check_species("CCO", ethanol("ethanol"), {"E": 1.0, "converged": False})
    # A fragment is used, flagged and kept out of the store
    result = {"E": 1.0, "converged": False}
    node._check_species("CCO", ethanol("CBH:CCO"), result)
    assert node._flags({"CCO": result})["CBH unconverged species"] == "CCO"
    assert node._storable("CCO", result) is False


def test_mopac_levels(node):
    node.model = MOPAC
    assert node._is_mopac(None)
    assert node._is_mopac({"level": "MOPAC:SQM@PM6-ORG"})
    assert not node._is_mopac({"level": "ORCA:DFT@r2SCAN-3c"})
    node.model = "ORCA:DFT@r2SCAN-3c"
    assert not node._is_mopac(None)


def test_below_cbh2_mopac_reports_its_own_value(node):
    """The fixture's best rung is CBH-1 (C + CCF -> CC + CF)."""
    node.model = MOPAC
    data = node.run_cbh(_parameters(node))
    assert data["CBH reaction"] == cbh_approach.MOPAC_OWN_LABEL
    assert data["CBH rung"] == 1
    # MOPAC's energy of the target is its enthalpy of formation
    assert data["DfH 298K"] == pytest.approx(ENERGIES["CCF"])
    # The CBH-1 value is kept for comparison
    assert data["CBH DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert "DfH uncertainty" not in data and "DfH 0K" not in data
    assert data["DfH uncertainty method"].startswith("none")


def test_other_levels_keep_cbh1(node):
    node.model = "ORCA:DFT@r2SCAN-3c"
    data = node.run_cbh(_parameters(node))
    assert data["CBH reaction"] == "CBH-1: C + CCF -> CC + CF"
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert "CBH DfH 298K" not in data


def test_no_rung_with_mopac_computes_the_target(node, monkeypatch):
    node.model = MOPAC
    system, configuration = node.get_system_configuration()
    computed, written = [], []

    def compute(P, can, smiles, system, configuration, embed=False):
        computed.append(can)
        return {"E": -123.4, "H298": -123.4}

    monkeypatch.setattr(node, "_compute_species", compute)
    monkeypatch.setattr(node, "_write_back", lambda P, c: written.append(list(c)))
    data = node._mopac_only_result(
        _parameters(node),
        "CCF",
        "CCF",
        system,
        configuration,
        SimpleNamespace(version="9.999"),
    )
    assert computed == ["CCF"] and written == [["CCF"]]
    assert data["DfH 298K"] == pytest.approx(-123.4)
    assert data["CBH reaction"] == cbh_approach.MOPAC_OWN_LABEL
    assert "CBH DfH 298K" not in data and "CBH rung" not in data
    assert data["CBH unconverged species"] == ""


def _reactions():
    """CBH-1 needs CF; CBH-2 needs CF and CCC; CBH-3 needs neither."""
    return {
        1: SimpleNamespace(coefficients={"T": -1, "C": -1, "CF": 1, "CC": 1}),
        2: SimpleNamespace(coefficients={"T": -1, "CF": 1, "CCC": 1}),
        3: SimpleNamespace(coefficients={"T": -1, "CCO": 1}),
    }


def test_an_unconverged_fragment_fails_its_rungs(node):
    energies = {
        "T": {"E": 1.0, "converged": False},  # the target is not a fragment
        "C": {"E": 1.0, "converged": True},
        "CF": {"E": 1.0, "converged": False},
        "CC": {"E": 1.0},
        "CCC": {"E": 1.0, "converged": False},
        "CCO": {"E": 1.0, "converged": True},
    }
    failed = node._failed_rungs("T", [1, 2, 3], _reactions(), energies)
    assert failed == {1: ["CF"], 2: ["CCC", "CF"]}
    assert node._failed_rungs_text(failed) == (
        "CBH-1: CF did not converge; CBH-2: CCC, CF did not converge"
    )
    # Converged in the extra steps: nothing fails
    energies["CF"]["converged"] = energies["CCC"]["converged"] = True
    assert node._failed_rungs("T", [1, 2, 3], _reactions(), energies) == {}


def test_a_failed_rung_falls_back_to_the_next_lower(node, monkeypatch):
    """The fixture's CBH-1 needs CF; with CBH-0 as well, CF failing leaves
    CBH-0."""
    node.model = "ORCA:DFT@r2SCAN-3c"
    real = node._energies_from_properties

    def energies(P, configuration, can):
        result = real(P, configuration, can)
        if can == "CF":
            result["converged"] = False
        return result

    monkeypatch.setattr(node, "_energies_from_properties", energies)
    # CBH-0 (CCF + 2 H2 -> 2 C + F) needs HF and H2 as well
    from tests.conftest import PROPERTY

    db = node._test_db
    current = db.system
    for smiles, E in (("F", -280.0), ("[H][H]", -10.0)):
        system = db.create_system(name=f"CBH:{smiles}", make_current=True)
        configuration = system.create_configuration(name="initial", make_current=True)
        configuration.from_smiles(smiles)
        configuration.properties.add(PROPERTY, "float", units="kJ/mol", noerror=True)
        configuration.properties.put(PROPERTY, E)
    db.system = current
    data = node.run_cbh(_parameters(node, **{"cbh rungs": "0, 1"}))
    assert data["CBH rung"] == 0
    assert data["CBH failed rungs"] == "CBH-1: CF did not converge"
    assert data["CBH unconverged species"] == "CF"


def test_with_mopac_a_failed_last_rung_gives_its_own_value(node, monkeypatch):
    node.model = MOPAC
    real = node._energies_from_properties

    def energies(P, configuration, can):
        result = real(P, configuration, can)
        if can == "CF":
            result["converged"] = False
        return result

    monkeypatch.setattr(node, "_energies_from_properties", energies)
    data = node.run_cbh(_parameters(node))
    assert data["CBH reaction"] == cbh_approach.MOPAC_OWN_LABEL
    assert data["DfH 298K"] == pytest.approx(ENERGIES["CCF"])
    assert data["CBH failed rungs"] == "CBH-1: CF did not converge"


def test_convergence_is_recorded_for_the_harvest(node, ethanol):  # noqa: F811
    node.model = MOPAC
    name = f"optimization converged#Thermochemistry#{MOPAC}"
    for result, P, expected in (
        ({"converged": True}, {}, "yes"),
        ({"converged": False}, {}, "no"),
        ({}, {"optimize fragments": "no"}, "not optimized"),
        ({"converged": None}, {"optimize fragments": "yes"}, "unknown"),
    ):
        configuration = ethanol("CBH:CCO")
        node._record_convergence(P, configuration, result)
        assert configuration.properties.get(name)[name]["value"] == expected
