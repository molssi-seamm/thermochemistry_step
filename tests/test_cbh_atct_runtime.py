# -*- coding: utf-8 -*-

"""ATcT at run time (P5d): the current version is downloaded when none is
installed, and covariance blocks are fetched by default, with quadrature when a
fetch fails. atct.anl.gov is never contacted: the download is the synthetic
fixture and the fetches are stand-ins."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from seamm_thermochemistry import atct as A  # noqa: E402

from tests.conftest import FIXTURE  # noqa: E402
from tests.test_cbh_existing import _parameters  # noqa: E402

QUADRATURE = (0.29**2 + 0.14**2 + 0.043**2) ** 0.5


@pytest.fixture()
def from_fixture(monkeypatch):
    """download_all from the synthetic fixture, counting the calls."""
    real = A.download_all
    calls = []

    def download(base=None, **kwargs):
        calls.append(base)
        return real(base, source=FIXTURE, **kwargs)

    monkeypatch.setattr(A, "download_all", download)
    return calls


def test_atct_is_downloaded_when_none_is_installed(node, tmp_path, from_fixture):
    empty = tmp_path / "no-ATcT"
    P = _parameters(node, **{"ATcT version": "latest", "ATcT directory": str(empty)})
    data = node.run_cbh(P)
    assert len(from_fixture) == 1
    assert data["ATcT version"] == "9.999"
    assert (empty / "9.999" / "atct_all.json").exists()
    # Installed now: not downloaded again
    node.run_cbh(P)
    assert len(from_fixture) == 1


def test_an_installed_version_is_never_downloaded(node, from_fixture):
    node.run_cbh(_parameters(node, **{"ATcT version": "latest"}))
    assert from_fixture == []


def test_a_specific_version_is_never_downloaded(node, tmp_path, from_fixture):
    P = _parameters(node, **{"ATcT directory": str(tmp_path / "no-ATcT")})
    with pytest.raises(RuntimeError, match="9.999 is not installed"):
        node.run_cbh(P)
    assert from_fixture == []


def test_a_failed_download_names_the_installer(node, tmp_path, monkeypatch):
    def refused(base=None, **kwargs):
        raise A.ATcTError("ATcT refused the request: HTTP 403.")

    monkeypatch.setattr(A, "download_all", refused)
    P = _parameters(
        node, **{"ATcT version": "latest", "ATcT directory": str(tmp_path / "none")}
    )
    with pytest.raises(RuntimeError) as e:
        node.run_cbh(P)
    assert "HTTP 403" in str(e.value)
    assert "seamm-thermochemistry-installer install" in str(e.value)


def test_covariances_are_fetched_by_default(node, monkeypatch):
    from thermochemistry_step.thermochemistry_parameters import (
        ThermochemistryParameters,
    )

    assert ThermochemistryParameters.parameters["covariances"]["default"] == (
        "fetch if needed"
    )
    fetched = []

    def fetch(self, ids, **kwargs):
        # A diagonal block: the same value as quadrature, but "fetched"
        u = {_id: self.entries[_id]["uncertainty_298K"] or 0.0 for _id in ids}
        block = {
            "ids": list(ids),
            "matrix": [[u[a] ** 2 if a == b else 0.0 for b in ids] for a in ids],
        }
        fetched.append(list(ids))
        self.covariance_cache["blocks"].append(block)
        return block

    monkeypatch.setattr(A.ATcT, "fetch_covariance", fetch)
    data = node.run_cbh(_parameters(node, covariances="fetch if needed"))
    assert len(fetched) == 1
    assert data["DfH uncertainty method"] == "ATcT covariance"
    assert data["DfH uncertainty"] == pytest.approx(QUADRATURE, abs=1e-6)


def test_a_failed_fetch_falls_back_to_quadrature(node, monkeypatch, capsys):
    attempts = []

    def offline(self, ids, **kwargs):
        attempts.append(list(ids))
        raise A.ATcTError("Cannot reach ATcT: no route to host")

    monkeypatch.setattr(A.ATcT, "fetch_covariance", offline)
    data = node.run_cbh(_parameters(node, covariances="fetch if needed"))
    assert len(attempts) == 1
    assert data["DfH uncertainty method"] == "quadrature (covariance fetch failed)"
    assert data["DfH uncertainty"] == pytest.approx(QUADRATURE, abs=1e-6)


def test_cached_only_never_fetches(node, monkeypatch):
    def never(self, ids, **kwargs):
        raise AssertionError("fetched with 'cached only'")

    monkeypatch.setattr(A.ATcT, "fetch_covariance", never)
    data = node.run_cbh(_parameters(node))  # the helper sets "cached only"
    assert data["DfH uncertainty method"] == "quadrature (no covariance cached)"
