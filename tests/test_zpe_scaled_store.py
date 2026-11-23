# -*- coding: utf-8 -*-

"""The species' properties carry the ZPE scaled as the CBH uses it (pm8-data
layout v2.7.6 §2c; science, thermo and P5, 2026-11-23):

* ``ZPE`` is the scaled ZPE, ``ZPE uncorrected`` the ZPE from the frequencies as
  computed, and H(0 K), H(298.15 K), E thermal and H thermal include the scaled
  ZPE; the scale factor and its source are stored with them;
* the readers that reuse them (a species already in the job, "existing results")
  unscale them, so the run's own factor is applied once: old (unscaled, no
  scale factor) and new (scaled) properties give the same enthalpy of formation.
"""

import pytest

from seamm_util import Q_

from tests.test_cbh_reuse import target_sets_fallback  # noqa: F401,E402
from tests.test_cbh_energy_level import split  # noqa: F401,E402
from .test_frequency_correction import (  # noqa: F401
    CH4_PM7,
    CH4_XYZ,
    _P,
    harmonic,
)

S = 0.98
RT = Q_(8.314462618 * 298.15, "J/mol").m_as("kJ/mol")


def _values(configuration, model):
    p = configuration.properties.get("*")
    return {
        k.split("#")[0]: v["value"]
        for k, v in p.items()
        if k.endswith(f"#Thermochemistry#{model}")
    }


@pytest.mark.parametrize("model", ["MOPAC:SQM@PM6-ORG", "MOPAC:SQM@PM7"])
def test_record_fields(harmonic, model):  # noqa: F811
    """Thermo's record checks (read_e1t_cbh.py), and science's: the stored H
    thermal is the H - E the CBH reaction uses, to 1e-6 kJ/mol."""
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7, model=model)
    corrected = model.endswith("PM7")
    node._cbh_P = {} if corrected else {"ZPE scale factor": str(S)}
    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    data["total_energy"] = -40.0
    configuration = node._working_configuration
    node._store_species_properties(configuration, data)
    v = _values(configuration, model)
    E = Q_(-40.0, "E_h").m_as("kJ/mol")
    s = 1.0 if corrected else S
    assert v["ZPE scale factor"] == s
    assert v["ZPE scale factor source"]
    # The H - E the CBH uses for this species: H thermal + (s - 1) ZPE
    used = data["H thermal"][0][0] + (s - 1.0) * data["ZPE"]
    assert v["H thermal(298.15 K)"] == pytest.approx(used, abs=1e-6)
    assert v["ZPE"] == pytest.approx(s * data["ZPE"], abs=1e-9)
    assert v["H(0 K)"] - E == pytest.approx(v["ZPE"], abs=1e-6)
    assert v["H(298.15 K)"] - E == pytest.approx(v["H thermal(298.15 K)"], abs=1e-6)
    assert v["H thermal(298.15 K)"] - v["E thermal(298.15 K)"] == pytest.approx(
        RT, abs=1e-3
    )
    if corrected:
        # PM7: corrected frequencies, no scaling; ZPE uncorrected from the raw ones
        assert v["ZPE uncorrected"] == pytest.approx(data["ZPE uncorrected"])
        assert v["ZPE uncorrected"] != pytest.approx(v["ZPE"])
        assert "frequency correction" in v
    else:
        assert v["ZPE"] / v["ZPE uncorrected"] == pytest.approx(S, abs=1e-12)


def test_unscale():
    from thermochemistry_step.cbh_approach import _unscale

    assert _unscale(9.8, 114.8, None) == (9.8, 114.8)  # old: stored unscaled
    ZPE, H = _unscale(9.8, 114.8, S)
    assert ZPE == pytest.approx(10.0) and H == pytest.approx(115.0)


# "Existing results (run nothing)" with harmonic corrections: per-species ZPEs
# that do not cancel in CBH-1 (CCF + C -> CF + CC), so double scaling would show
ZPES = {"CCF": 160.0, "C": 110.0, "CF": 100.0, "CC": 150.0}
H_THERMAL = {can: z + 12.0 for can, z in ZPES.items()}  # unscaled H - E
MODEL = "Fake:M"


def _existing(node, monkeypatch, style):
    from tests.conftest import PROPERTY
    from tests.test_cbh_existing import _parameters

    for system in node._test_db.get_systems():
        can = system.name[4:]
        configuration = system.configuration
        p = configuration.properties
        E = p.get(PROPERTY)[PROPERTY]["value"]
        z, h = ZPES[can], H_THERMAL[can]
        if style == "new":  # as stored since 2026.11.23
            shift = (S - 1.0) * z
            values = {"ZPE": z + shift, "H(298.15 K)": E + h + shift}
            values["ZPE uncorrected"] = z
            values["ZPE scale factor"] = S
        else:  # before: unscaled, no scale factor (jobs 5955, 5962)
            values = {"ZPE": z, "H(298.15 K)": E + h}
        for base, value in values.items():
            name = f"{base}#Thermochemistry#{MODEL}"
            units = "" if base == "ZPE scale factor" else "kJ/mol"
            p.add(name, "float", units=units, noerror=True)
            p.put(name, value)
    return _parameters(
        node,
        **{
            "thermal corrections": "harmonic",
            "energy property": "total energy#Test#*",
            "low-level energy property": "total energy#Test#*",
            "ZPE scale factor": str(S),
        },
    )


@pytest.mark.parametrize("style", ["old", "new"])
def test_existing_results_scale_once(node, monkeypatch, style):
    pytest.importorskip("rdkit")
    P = _existing(node, monkeypatch, style)
    data = node.run_cbh(P)
    # DrH = sum of nu * (E + h + (S - 1) z) over CF + CC - CCF - C
    nu = {"CF": 1, "CC": 1, "CCF": -1, "C": -1}
    thermal = sum(n * (H_THERMAL[c] + (S - 1.0) * ZPES[c]) for c, n in nu.items())
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0 - thermal, abs=1e-6)
    # Not the doubly scaled value
    assert abs(sum(n * (S - 1.0) ** 2 * ZPES[c] for c, n in nu.items())) > 1e-3


# Species already computed in the job (an earlier Loop iteration)
def test_in_job_species_scale_once(target_sets_fallback):  # noqa: F811
    """Methane is in the job with old (unscaled) properties; the first run reuses
    it and computes the rest. Then every species is stored the new (scaled) way,
    and the second run reuses them all: the same enthalpy of formation, so
    nothing is scaled twice."""
    pytest.importorskip("rdkit")
    from tests.test_cbh_energy_level import LOW, FakeEvaluator, _computed
    from tests.test_cbh_reuse import _computed_in_job, _put

    node = target_sets_fallback
    configuration = node._find_fragment(node._test_db, "C")
    _computed_in_job(configuration, "C", "fp:low", [300.0, 1500.0])
    s = f"#Thermochemistry#{LOW}"
    p = configuration.properties
    E = p.get("H(298.15 K)" + s)["H(298.15 K)" + s]["value"] - 15.0
    z, h = 100.0, 115.0
    _put(configuration, "ZPE" + s, "float", z, "kJ/mol")
    _put(configuration, "H(298.15 K)" + s, "float", E + h, "kJ/mol")
    P = _computed(node, **{"ZPE scale factor": str(S)})
    FakeEvaluator.calls = []
    first = node.run_cbh(P)["DfH 298K"]
    assert "C" not in node._computed  # reused, not computed
    computed = set(node._computed) - {"CCF"}
    assert computed
    # The computed species as the step now stores them: scaled, with the factor
    # (the stub computed them with ZPE 10 and H - E 15)
    for can in computed:
        other = node._find_fragment(node._test_db, can)
        _computed_in_job(other, can, "fp:low", [300.0, 1500.0])
        E_can = other.properties.get("H(298.15 K)" + s)["H(298.15 K)" + s]["value"]
        E_can -= 15.0
        _put(other, "ZPE" + s, "float", S * 10.0, "kJ/mol")
        _put(other, "H(298.15 K)" + s, "float", E_can + 15.0 + (S - 1) * 10.0, "kJ/mol")
        _put(other, "ZPE scale factor" + s, "float", S, "")
    # Methane the new way
    _put(configuration, "ZPE" + s, "float", S * z, "kJ/mol")
    _put(configuration, "H(298.15 K)" + s, "float", E + h + (S - 1) * z, "kJ/mol")
    _put(configuration, "ZPE scale factor" + s, "float", S, "")
    FakeEvaluator.calls = []
    node._computed.clear()
    second = node.run_cbh(P)["DfH 298K"]
    # Every fragment reused (the target is computed again, as always)
    assert set(node._computed) <= {"CCF"}, node._computed
    assert second == pytest.approx(first, abs=1e-6)
