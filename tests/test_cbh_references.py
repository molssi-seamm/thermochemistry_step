# -*- coding: utf-8 -*-

"""The supplementary reference library in the CBH approach: a reference that
ATcT lacks comes from the library, its uncertainty is added in quadrature, and
the result says where its references came from."""

import pytest

pytest.importorskip("rdkit")
references = pytest.importorskip("seamm_thermochemistry.references")

from seamm_thermochemistry import atct as A  # noqa: E402

from tests.test_cbh_existing import _parameters, radical_node  # noqa: E402, F401


@pytest.fixture()
def no_cf(monkeypatch):
    """ATcT without fluoromethane, so CBH-1 of fluoroethane needs the library."""
    real = A.ATcT.lookup

    def lookup(self, smiles=None, **kwargs):
        if smiles is not None and A.canonical_smiles(smiles) == "CF":
            raise A.MissingEntry(f"No ATcT entry for {smiles}")
        return real(self, smiles, **kwargs)

    monkeypatch.setattr(A.ATcT, "lookup", lookup)


def _library(node, value=-234.0, u=1.0, offset=None):
    where = node._test_base.parent
    lib = references.ReferenceLibrary(
        where / "references-official.db", create=True, layer="official"
    )
    lib.add_source("WebBook-SRD69", "measured", 50)
    sid = lib.add_species("CF", formula="CH3F")
    lib.add_measurement(
        sid,
        "WebBook-SRD69",
        dfh_298k=value,
        uncertainty=u,
        unc_level="2 sigma",
        citation="Someone et al. 1990",
        year=1990,
    )
    lib.choose(sid, "WebBook-SRD69")
    if offset is not None:
        lib.set_offset("WebBook-SRD69", offset, 9, "9.999")
    lib.close()


def _printed(caplog):
    return "\n".join(r.getMessage() for r in caplog.records)


# Fluoroethane is itself in ATcT, so without a usable rung its own value is
# the result; what the reactions report shows where the references were sought


def test_atct_only_without_the_library(node, no_cf, caplog):
    node.run_cbh(_parameters(node, **{"supplementary references": "none"}))
    assert "not usable: no ATcT entry for CF" in _printed(caplog)
    assert "reference library" not in _printed(caplog)


def test_a_missing_library_falls_back_to_atct(node, no_cf, caplog):
    node.run_cbh(_parameters(node, **{"supplementary references": "compilations"}))
    assert "No supplementary reference library is installed" in _printed(caplog)
    assert "not usable: no ATcT entry for CF" in _printed(caplog)


def test_reference_from_the_library(node, no_cf):
    _library(node, offset=1.16)
    data = node.run_cbh(
        _parameters(node, **{"supplementary references": "compilations"})
    )
    assert data["CBH rung"] == 1
    # CF -234.0 (library) + CC -83.91 - C -74.513 (ATcT), DrH 2.0
    assert data["DfH 298K"] == pytest.approx(-234.0 - 83.91 + 74.513 - 2.0)
    u_lib = (1.0**2 + 1.16**2) ** 0.5
    assert data["DfH uncertainty"] == pytest.approx(
        (0.14**2 + 0.043**2 + u_lib**2) ** 0.5, abs=1e-6
    )
    assert data["DfH uncertainty method"].endswith("+ reference library (quadrature)")
    assert data["CBH references"] == "mixed"
    used = data["references"]["library values"]
    assert used["CF"]["source"] == "WebBook-SRD69"
    assert used["CF"]["citation"] == "Someone et al. 1990"
    assert data["references"]["ATcT"] == "9.999"


def test_the_library_is_not_used_when_atct_has_the_species(node):
    _library(node, value=-100.0)
    data = node.run_cbh(
        _parameters(node, **{"supplementary references": "compilations"})
    )
    assert data["CBH references"] == "ATcT-only"
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)


def test_no_value_anywhere(node, no_cf, caplog):
    lib = references.ReferenceLibrary(
        node._test_base.parent / "references-official.db",
        create=True,
        layer="official",
    )
    lib.add_source("JANAF-1998", "evaluated", 10)
    lib.close()
    node.run_cbh(_parameters(node, **{"supplementary references": "compilations"}))
    assert "not usable: no reference value for CF" in _printed(caplog)


def test_parameters():
    from thermochemistry_step.thermochemistry_parameters import (
        ThermochemistryParameters,
    )

    p = ThermochemistryParameters.parameters
    assert p["supplementary references"]["default"] == "compilations and computed"
    assert p["use my own references"]["default"] == "no"


def test_the_error_names_the_remedies(radical_node, monkeypatch):  # noqa: F811
    from tests.test_cbh_existing import _parameters as parameters

    node = radical_node
    monkeypatch.setattr(node._test_target, "to_smiles", lambda **k: "[CH2]CC")
    real = A.ATcT.lookup

    def lookup(self, smiles=None, **kwargs):
        if smiles is not None and A.canonical_smiles(smiles) == "[CH2]C":
            raise A.MissingEntry(f"No ATcT entry for {smiles}")
        return real(self, smiles, **kwargs)

    monkeypatch.setattr(A.ATcT, "lookup", lookup)
    lib = references.ReferenceLibrary(
        node._test_base.parent / "references-official.db",
        create=True,
        layer="official",
    )
    lib.add_source("JANAF-1998", "evaluated", 10)
    lib.close()
    P = parameters(node, **{"supplementary references": "compilations"})
    with pytest.raises(RuntimeError, match="or the reference library") as e:
        node.run_cbh(P)
    assert "Use my own references" in str(e.value)
