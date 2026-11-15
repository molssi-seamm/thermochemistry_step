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
