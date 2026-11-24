# -*- coding: utf-8 -*-

"""The thermal corrections at 298.15 K and the ZPE scale factor as properties
(pm8-data layout §2c, science 2026-11-23): ``ZPE scale factor``, ``ZPE scale factor
source``, ``E thermal(298.15 K)`` and ``H thermal(298.15 K)``, keyed
``#Thermochemistry#{model}``, in kJ/mol where they have units. E thermal = U(T) - E
and H thermal = H(T) - E, both including the ZPE."""

import csv
from pathlib import Path

import pytest

from seamm_util import Q_

from .test_frequency_correction import CH4_PM7, CH4_XYZ, _P, harmonic  # noqa: F401

CSV = Path(__file__).resolve().parents[1] / "thermochemistry_step" / "data"
NAMES = (
    "ZPE scale factor",
    "ZPE scale factor source",
    "E thermal(298.15 K)",
    "H thermal(298.15 K)",
)


def test_scalars_at_298_and_1_bar(harmonic):  # noqa: F811
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    assert data["E thermal(298.15 K)"] == data["E thermal"][0][0]
    assert data["H thermal(298.15 K)"] == data["H thermal"][0][0]
    # Both include the ZPE, and H - U = RT for an ideal gas
    assert data["E thermal(298.15 K)"] > data["ZPE"]
    RT = Q_(8.314462618 * 298.15, "J/mol").m_as("kJ/mol")
    assert data["H thermal(298.15 K)"] - data["E thermal(298.15 K)"] == pytest.approx(
        RT, abs=1e-3
    )


def test_picked_from_several_temperatures(harmonic):  # noqa: F811
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    data = node.harmonic_thermochemistry(
        _P(), [200.0, 298.15, 400.0], "K", [2.0, 1.0], "bar"
    )
    assert data["E thermal(298.15 K)"] == data["E thermal"][1][1]
    assert data["H thermal(298.15 K)"] == data["H thermal"][1][1]


@pytest.mark.parametrize("Ts, Ps", [([300.0], [1.0]), ([298.15], [2.0])])
def test_absent_without_298_and_1_bar(harmonic, Ts, Ps):  # noqa: F811
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    data = node.harmonic_thermochemistry(_P(), Ts, "K", Ps, "bar")
    assert "E thermal(298.15 K)" not in data
    assert "H thermal(298.15 K)" not in data


def test_cbh_species_carry_the_corrections(harmonic):  # noqa: F811
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    data["total_energy"] = -40.0
    configuration = node._working_configuration
    node._store_species_properties(configuration, data)
    values = configuration.properties.get("*")
    model = node.model
    E = values[f"E thermal(298.15 K)#Thermochemistry#{model}"]["value"]
    H = values[f"H thermal(298.15 K)#Thermochemistry#{model}"]["value"]
    assert E == pytest.approx(data["E thermal"][0][0])
    assert H == pytest.approx(data["H thermal"][0][0])
    # Consistent with the stored absolute enthalpy
    H298 = values[f"H(298.15 K)#Thermochemistry#{model}"]["value"]
    assert H298 == pytest.approx(Q_(-40.0, "E_h").m_as("kJ/mol") + H)
    for name in ("E thermal(298.15 K)", "H thermal(298.15 K)"):
        assert configuration.properties.units(f"{name}#Thermochemistry#{model}") == (
            "kJ/mol"
        )


def test_definitions():
    from thermochemistry_step.metadata import metadata

    results = metadata["results"]
    rows = {
        row["Property"]: row
        for row in csv.DictReader(open(CSV / "properties.csv", encoding="utf-8"))
    }
    for name in NAMES:
        assert results[name]["property"] == f"{name}#Thermochemistry#{{model}}"
        row = rows[f"{name}#Thermochemistry#{{model}}"]
        if "thermal" in name:
            assert results[name]["units"] == "kJ/mol" and row["Units"] == "kJ/mol"
            assert row["Type"] == "float"
    assert rows["ZPE scale factor source#Thermochemistry#{model}"]["Type"] == "str"


# <S**2> of open-shell species from the warm engine (item 1 of 8e49882f)
def _species(harmonic):  # noqa: F811
    node = harmonic(["C", "H", "H", "H", "H"], CH4_XYZ, CH4_PM7)
    data = node.harmonic_thermochemistry(_P(), [298.15], "K", [1.0], "bar")
    data["total_energy"] = -40.0
    return node, data


def test_spin_stored_for_an_open_shell(harmonic):  # noqa: F811
    node, data = _species(harmonic)
    node._spin = (0.7556, 0.75)
    configuration = node._working_configuration
    node._store_species_properties(configuration, data)
    values = configuration.properties.get("*")
    assert values[f"S^2#ORCA#{node.model}"]["value"] == pytest.approx(0.7556)
    assert values[f"S^2 ideal#ORCA#{node.model}"]["value"] == pytest.approx(0.75)
    assert node._spin is None  # taken once: not given to the next species


@pytest.mark.parametrize("spin", [None, (None, None)])
def test_no_spin_for_a_closed_shell_or_without_it(harmonic, spin):  # noqa: F811
    node, data = _species(harmonic)
    node._spin = spin
    configuration = node._working_configuration
    node._store_species_properties(configuration, data)
    names = configuration.properties.get("*").keys()
    assert not any(name.startswith("S^2") for name in names)


def test_the_driver_spin_reaches_the_step(node, monkeypatch):
    """The Hessian mixin copies the driver's final <S**2>, like its fingerprints."""

    class Driver:
        fingerprints = None
        last_optimization = None
        spin = (2.01, 2.0)

        def optimize_and_energy(self, *args, **kwargs):
            return -1.0

    monkeypatch.setattr(node, "_hessian_driver", lambda: Driver())
    node.optimize_and_energy(None)
    assert node._spin == (2.01, 2.0)
