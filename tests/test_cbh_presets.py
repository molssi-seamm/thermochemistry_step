# -*- coding: utf-8 -*-

"""The presets of the CBH approach (quick, default, custom), the ZPE scale
factor, the spread between the rungs and the warning of the quick preset."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

import seamm_exec  # noqa: E402

from thermochemistry_step import cbh_approach  # noqa: E402
from tests.conftest import ENERGIES  # noqa: E402
from tests.test_cbh_existing import _parameters  # noqa: E402

R2SCAN = "ORCA:CMP@R2SCAN-3C"
REVDSD = "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP"
ZPE = 10.0


class FakeEvaluator:
    """The energy level's energies: ENERGIES (by system name)."""

    calls = []

    def __init__(self, node, mc, **kwargs):
        FakeEvaluator.level = mc["level"]
        self.submitted = []

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def submit(self, configuration, key):
        self.submitted.append((key, configuration))

    def results(self):
        names = []
        for key, configuration in self.submitted:
            name = configuration.system.name
            names.append(name)
            can = name[4:] if name.startswith("CBH:") else "CCF"
            yield seamm_exec.EvaluatorResult(
                key=key, ok=True, energy=ENERGIES.get(can, -100.0)
            )
        FakeEvaluator.calls.append(names)


@pytest.fixture()
def preset(node, monkeypatch):
    """The node with the presets' model chemistries resolved by name, its
    species computed (stubbed: the geometry level's energy is ENERGIES + 50,
    ZPE 10, H - E 15) and the printed text collected."""
    FakeEvaluator.calls = []
    node.model = R2SCAN
    monkeypatch.setattr(seamm_exec, "Evaluator", FakeEvaluator)
    monkeypatch.setattr(
        node, "_resolve_model_chemistry", lambda text: {"level": text, "step": "X"}
    )

    def compute(P, can, generated, system=None, configuration=None, embed=False):
        # As the real one does, a fragment gets its system CBH:<SMILES>
        db = node._test_db
        if node._find_fragment(db, can) is None:
            new = db.create_system(name=f"CBH:{can}", make_current=False)
            new.create_configuration(name="initial", make_current=True).from_smiles(
                generated
            )
        E = ENERGIES.get(can, -100.0) + 50.0
        return {"E": E, "ZPE": ZPE, "H0": E + ZPE, "H298": E + 15.0}

    monkeypatch.setattr(node, "_compute_species", compute)
    printed = []
    monkeypatch.setattr(
        cbh_approach.printer,
        "important",
        lambda text, *a, **k: printed.append(str(text)),
    )
    node._printed = printed
    return node


def _preset(node, name, **changes):  # noqa: F811
    return _parameters(node, preset=name, **changes)


def test_apply_preset(node):
    P = _parameters(node, preset="default")
    Q = node.apply_preset(P)
    assert Q["cbh rungs"] == "0, 1, 2" and Q["thermal corrections"] == "harmonic"
    assert Q["cbh energies"] == "finite-difference harmonic"
    assert Q["energy model chemistry"] == REVDSD
    assert node.apply_preset({**P, "preset": "quick"})["energy model chemistry"] == (
        "the geometry level"
    )
    assert node.apply_preset({**P, "include CBH-3": "yes"})["cbh rungs"] == "0, 1, 2, 3"
    # custom leaves everything alone
    custom = _parameters(node, preset="custom")
    assert node.apply_preset(custom) is custom


def test_preset_model_chemistry(preset):
    assert preset.preset_model_chemistry({"preset": "default"})["level"] == R2SCAN
    assert preset.preset_model_chemistry({"preset": "custom"}) is None


def test_zpe_scale(node):
    node.model = R2SCAN
    factor, source, flag = node._zpe_scale({"ZPE scale factor": "from the table"})
    assert factor == 0.9899 and "W4-2006" in source and "extrapolation" in flag
    assert node._zpe_scale({"ZPE scale factor": "0.97"}) == (0.97, "as given", "")
    node.model = "Test:PM7"
    factor, source, flag = node._zpe_scale({})
    assert factor == 1.0 and "no ZPE scale factor" in flag
    with pytest.raises(ValueError, match="must be"):
        node._zpe_scale({"ZPE scale factor": "lots"})


def test_default_preset(preset):
    data = preset.run_cbh(_preset(preset, "default"))
    # revDSD energies (ENERGIES) at the highest usable rung, CBH-1 here; the ZPE
    # scaling and the thermal part cancel in this reaction (2 molecules each side)
    assert data["CBH rung"] == 1
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert data["CBH energy level"] == REVDSD and data["CBH geometry level"] == R2SCAN
    assert FakeEvaluator.level == REVDSD
    assert data["ZPE scale factor"] == 0.9899
    assert "extrapolation" in data["ZPE scale factor source"]
    assert any("The 'default' preset" in t for t in preset._printed)
    assert not any("Use the 'default' preset" in t for t in preset._printed)


def test_quick_preset_warns_below_cbh2(preset):
    data = preset.run_cbh(_preset(preset, "quick"))
    assert FakeEvaluator.calls == []  # the geometry level's own energies
    assert data["CBH energy level"] == R2SCAN
    assert data["CBH rung"] == 1
    assert any("Use the 'default' preset" in t for t in preset._printed)


def test_zpe_scale_is_applied(preset):
    """The factor reaches the result; in CCF + C -> CF + CC (two molecules each
    side, equal ZPEs) the shifts cancel."""
    custom = _preset(
        preset,
        "custom",
        **{
            "cbh energies": "finite-difference harmonic",
            "thermal corrections": "harmonic",
            "energy model chemistry": "the geometry level",
            "ZPE scale factor": "0.5",
        },
    )
    data = preset.run_cbh(custom)
    assert data["ZPE scale factor"] == 0.5
    unscaled = preset.run_cbh({**custom, "ZPE scale factor": "1.0"})
    assert data["DfH 0K"] == pytest.approx(unscaled["DfH 0K"])


def test_rung_spread():
    results = {
        0: {"DfH 298.15 K": -10.0},
        1: {"DfH 298.15 K": -100.0},
        2: {"DfH 298.15 K": -101.5},
        3: {"DfH 298.15 K": -100.5},
    }
    spread, rungs = cbh_approach.rung_spread(results, [0, 1, 2, 3])
    assert rungs == [1, 2, 3] and spread == pytest.approx(1.5)  # CBH-0 left out
    assert cbh_approach.rung_spread(results, [0, 1]) == (None, [1])


# -- the presets as composite recipes (composite_step) ------------------------------


@pytest.mark.parametrize("name", ["quick", "default"])
def test_presets_are_their_recipes(name):
    """The spellings in PRESETS (the fallback without composite_step) are the
    recipes' levels, and the recipes' ZPE factor is the table's."""
    pytest.importorskip("composite_step")
    from seamm_thermochemistry import zpe_scale_factor

    recipe = cbh_approach.recipe_levels(cbh_approach.PRESETS[name]["recipe"])
    assert recipe["geometry"] == cbh_approach.PRESETS[name]["geometry"]
    assert recipe["energy"] == cbh_approach.PRESETS[name]["energy"]
    assert recipe["zpe"][0] == zpe_scale_factor(recipe["geometry"]).factor


def test_default_preset_records_its_recipe(preset):
    pytest.importorskip("composite_step")
    P = _preset(preset, "default")
    preset.preset_model_chemistry(P)
    data = preset.run_cbh(P)
    assert data["CBH composite recipe"] == "SEAMM:CMP@revDSD-r2SCAN-3c-HO"
    assert data["CBH composite recipe version"] == 1
    assert data["CBH composite recipe fingerprint"].startswith("sha256:")
    # The levels, which key the species' cached properties and the store, are
    # the ordinary ones
    assert data["CBH geometry level"] == R2SCAN and preset.model == R2SCAN
    assert data["CBH energy level"] == REVDSD
    assert "extrapolation" in data["ZPE scale factor source"]  # the table's


def test_the_recipe_path_reproduces_the_preset_path(preset, monkeypatch):
    """With and without composite_step the preset gives the same levels and the
    same numbers: species cached or stored under the old presets are hit."""
    pytest.importorskip("composite_step")
    P = _preset(preset, "default")
    preset.preset_model_chemistry(P)
    with_recipe = preset.run_cbh(P)
    import builtins

    real_import = builtins.__import__

    def no_composite(name, *args, **kwargs):
        if name.startswith("composite_step"):
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_composite)
    preset.preset_model_chemistry(P)
    assert preset._composite is None
    without = preset.run_cbh(P)
    recipe_keys = {k for k in with_recipe if k.startswith("CBH composite recipe")}
    assert {k: v for k, v in with_recipe.items() if k not in recipe_keys} == without


def test_composite_as_the_flowcharts_model_chemistry(node, monkeypatch):
    pytest.importorskip("composite_step")
    import seamm

    monkeypatch.setattr(
        node, "_resolve_model_chemistry", lambda text: {"level": text, "step": "X"}
    )
    monkeypatch.setattr(seamm, "flowchart_variables", object())  # a flowchart
    variables = {
        "_model_chemistry": {"level": "SEAMM:CMP@revDSD-r2SCAN-3c-HO", "step": "C"}
    }
    monkeypatch.setattr(node, "variable_exists", lambda name: name in variables)
    monkeypatch.setattr(node, "get_variable", lambda name: variables[name])
    mc = node.preset_model_chemistry({"preset": "custom"})
    assert mc["level"] == R2SCAN
    node.model = mc["level"]
    energy = node._energy_level({"energy model chemistry": "the geometry level"})
    assert energy["level"] == REVDSD
    assert node._zpe_scale({"ZPE scale factor": "from the table"})[0] == 0.9899


def test_composite_as_the_energy_level(node, monkeypatch):
    pytest.importorskip("composite_step")
    monkeypatch.setattr(
        node, "_resolve_model_chemistry", lambda text: {"level": text, "step": "X"}
    )
    node._composite = None
    node.model = "ORCA:DFT@PBE0-D4/def2-TZVP"
    printed = []
    monkeypatch.setattr(
        cbh_approach.printer, "important", lambda t, *a, **k: printed.append(str(t))
    )
    energy = node._energy_level(
        {"energy model chemistry": "SEAMM:CMP@revDSD-r2SCAN-3c-HO"}
    )
    assert energy["level"] == REVDSD
    assert any("defined with geometries at" in " ".join(t.split()) for t in printed)


def test_multi_component_recipes_are_refused(tmp_path):
    pytest.importorskip("composite_step")
    (tmp_path / "recipes").mkdir()
    (tmp_path / "recipes" / "toy.yaml").write_text(
        "recipe: toy\nversion: 1\ngeometry: {level: 'ORCA:CMP@R2SCAN-3C'}\n"
        "components:\n  a: {level: 'ORCA:QC@MP2/def2-TZVP'}\n"
        "  b: {level: 'ORCA:QC@MP2/def2-SVP'}\ncombine: |\n  E0 = a - b\n"
    )
    with pytest.raises(RuntimeError, match="combines several calculations"):
        cbh_approach.recipe_levels("SEAMM:CMP@toy", data_path=[tmp_path])
