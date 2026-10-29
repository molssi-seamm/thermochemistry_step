# -*- coding: utf-8 -*-

"""The CBH approach in its "existing results (run nothing)" mode: the species'
energies are read from properties stored on the systems named CBH:<SMILES>,
and combined with ATcT references from a synthetic ATcT cache in tests/data
(never the real tables)."""

from pathlib import Path

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

import seamm  # noqa: E402
from molsystem import SystemDB  # noqa: E402
from seamm_thermochemistry import atct  # noqa: E402

import thermochemistry_step  # noqa: E402

FIXTURE = Path(__file__).parent / "data" / "atct_fixture.json"
PROPERTY = "total energy#Test#PM7"

# Synthetic energies, kJ/mol, taken as enthalpies at 298.15 K (no thermal
# corrections). Fluoroethane CBH-1: CCF + C -> CF + CC, so
# DrH = E(CF) + E(CC) - E(CCF) - E(C) = -230 - 80 + 262 + 50 = 2.0
ENERGIES = {"CCF": -262.0, "C": -50.0, "CF": -230.0, "CC": -80.0}


@pytest.fixture()
def atct_base(tmp_path):
    base = tmp_path / "ATcT"
    atct.download_all(base, source=FIXTURE)
    return base


@pytest.fixture()
def node(tmp_path, atct_base, monkeypatch):
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


def _parameters(node, **changes):
    P = {
        "cbh rungs": "1",
        "protect": "default",
        "cbh energies": "existing results (run nothing)",
        "thermal corrections": "none (the energy is the enthalpy, e.g. MOPAC)",
        "energy property": "total energy#*#*",
        "ZPE property": "ZPE#*#*",
        "enthalpy property": "H(298.15 K)#*#*",
        "low-level energy property": "",
        "hessian route": "automatic",
        "initial conformer": "lowest of the embeddings (MMFF)",
        "conformer embeddings": 20,
        "target conformer": "as given",
        "fragment store": "ignore",
        "fragment store file": "default",
        "reference values": "ATcT",
        "ladder": "",
        "ATcT uncertainty threshold": 0.5,
        "ATcT version": "9.999",
        "ATcT directory": str(node._test_base),
        "covariances": "cached only",
        "compare with ATcT": "yes",
    }
    P.update(changes)
    return P


def test_existing_results(node):
    data = node.run_cbh(_parameters(node))
    assert data["CBH rung"] == 1
    assert data["CBH reaction"] == "CBH-1: C + CCF -> CC + CF"
    # ATcT 298 K references of the fixture: CF -235.55 + CC -83.91 - C -74.513
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    assert "DfH 0K" not in data  # no 0 K value without thermal corrections
    assert data["DfH uncertainty"] == pytest.approx(
        (0.29**2 + 0.14**2 + 0.043**2) ** 0.5, abs=1e-6
    )
    assert data["ATcT DfH 298K"] == pytest.approx(-272.11)
    assert data["ATcT version"] == "9.999"
    table = seamm.Table(node._test_db, "CBH")
    assert table.n_rows == 1
    assert data["CBH results"]["1"]["uncertainty method"].startswith("quadrature")


# The n-propyl radical, CBH-1: [CH2]CC + C -> [CH2]C + CC (the ethyl
# radical's own CBH-1 is trivial), DrH = E([CH2]C) + E(CC) - E([CH2]CC) - E(C)
# = 130 - 80 - 110 + 50 = -10
RADICAL_ENERGIES = {"[CH2]CC": 110.0, "C": -50.0, "[CH2]C": 130.0, "CC": -80.0}


@pytest.fixture()
def radical_node(tmp_path, atct_base, monkeypatch):
    db = SystemDB(filename=f"file:{tmp_path / 'radical.db'}")
    target = None
    for smiles, energy in RADICAL_ENERGIES.items():
        system = db.create_system(name=f"CBH:{smiles}", make_current=True)
        configuration = system.create_configuration(name="initial", make_current=True)
        configuration.from_smiles(smiles)
        configuration.properties.add(PROPERTY, "float", units="kJ/mol", noerror=True)
        configuration.properties.put(PROPERTY, energy)
        if smiles == "[CH2]CC":
            target = (system, configuration)
    db.system = target[0]
    # from_smiles sets the multiplicity of a radical correctly ...
    assert target[1].spin_multiplicity == 2 and target[1].charge == 0

    node = thermochemistry_step.Thermochemistry()
    variables = {"_system_db": db}
    monkeypatch.setattr(node, "get_variable", lambda name: variables[name])
    monkeypatch.setattr(node, "get_system_configuration", lambda *a, **k: target)
    monkeypatch.setattr(node, "get_table", lambda name: seamm.Table.create(db, name))
    monkeypatch.setattr(node, "_cite_cbh", lambda tables: None)
    node._test_db = db
    node._test_base = atct_base
    node._test_target = target[1]
    yield node
    db.close()


def test_radical_target(radical_node, monkeypatch):
    node = radical_node
    # ... but to_smiles of a radical gives an anion (molsystem, 2026-10-28),
    # so stand in for the fixed conversion here
    monkeypatch.setattr(node._test_target, "to_smiles", lambda **k: "[CH2]CC")
    data = node.run_cbh(_parameters(node))
    from seamm_thermochemistry import cbh

    reaction = cbh.cbh_reactions("[CH2]CC", rungs=[1])[1]
    assert dict(reaction.references) == {"[CH2]C": 1, "CC": 1, "C": -1}
    assert data["CBH reaction"] == str(reaction)
    # 298 K references: [CH2]C 120.75 + CC -83.91 - C -74.513 = 111.353
    assert data["DfH 298K"] == pytest.approx(111.353 + 10.0, abs=1e-6)
    assert "ATcT DfH 298K" not in data  # the propyl radical is not in the fixture


def test_radical_smiles_guard(radical_node, monkeypatch):
    """A SMILES that disagrees with the configuration's charge and spin (as
    molsystem before 2026.10.28 gave for a radical, [CH2-]CC) must stop the
    step rather than build the reactions of the anion."""
    monkeypatch.setattr(radical_node._test_target, "to_smiles", lambda **k: "[CH2-]CC")
    with pytest.raises(RuntimeError, match=r"implies charge -1 and 0 unpaired"):
        radical_node.run_cbh(_parameters(radical_node))


def test_radical_target_from_structure(radical_node):
    """With molsystem's radical-aware SMILES (2026.10.28), the target's SMILES
    comes straight from the structure."""
    smiles = radical_node._test_target.to_smiles(canonical=True, isomeric=True)
    if "-" in smiles:
        pytest.skip(f"molsystem gives {smiles} for the n-propyl radical")
    data = radical_node.run_cbh(_parameters(radical_node))
    assert data["DfH 298K"] == pytest.approx(111.353 + 10.0, abs=1e-6)


def test_fragment_charge_and_multiplicity():
    from thermochemistry_step.cbh_approach import cbh_charge, cbh_multiplicity

    assert (cbh_charge("[CH3]"), cbh_multiplicity("[CH3]")) == (0, 2)
    assert (cbh_charge("[NH3+]C"), cbh_multiplicity("[NH3+]C")) == (1, 1)
    assert (cbh_charge("CC(=O)[O-]"), cbh_multiplicity("CC(=O)[O-]")) == (-1, 1)
    assert (cbh_charge("C[N+](=O)[O-]"), cbh_multiplicity("C[N+](=O)[O-]")) == (0, 1)
    assert (cbh_charge("[CH2]"), cbh_multiplicity("[CH2]")) == (0, 3)


def test_missing_fragment(node):
    with pytest.raises(RuntimeError, match="No property matching"):
        node.run_cbh(_parameters(node, **{"energy property": "nothing#*#*"}))


def test_fragment_from_store(node, monkeypatch, tmp_path):
    """A fragment with no system is served from the fragment store at the
    step's model chemistry; an unverified record is reported as such."""
    from seamm_thermochemistry import FragmentRecord, FragmentStore

    path = tmp_path / "fragments.db"
    with FragmentStore(path, create=True) as store:
        store.put(FragmentRecord.make("CF", E=-230.0, model_chemistry="Test:PM7"))
    find = node._find_fragment
    monkeypatch.setattr(
        node,
        "_find_fragment",
        lambda db, can: None if can == "CF" else find(db, can),
    )
    node.model = "Test:PM7"
    P = _parameters(node, **{"fragment store": "use if available"})
    P["fragment store file"] = str(path)
    data = node.run_cbh(P)
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0, abs=1e-6)
    # At another model chemistry the store has nothing
    node.model = "Other:PM6"
    with pytest.raises(RuntimeError, match="nothing in the fragment store"):
        node.run_cbh(P)


def _ladder(tmp_path, atct_base, shift):
    """A ladder whose CH3F and C2H6 values are ATcT's shifted by `shift`."""
    from seamm_thermochemistry import ATcT, build_ladder, FragmentRecord, FragmentStore

    tables = ATcT(base=atct_base)
    store_path = tmp_path / "fragments.db"
    level = "High:F12"
    with FragmentStore(store_path, create=True) as store:
        for smiles in ("[H][H]", "C", "F", "CC", "CF"):
            h0 = tables.lookup(smiles)["dfH0_0K"] + (
                shift if smiles in ("CC", "CF") else 0
            )
            store.put(
                FragmentRecord.make(
                    smiles, E=h0 - 10.0, ZPE=10.0, H=h0 + 5.0, model_chemistry=level
                )
            )
    ladder = build_ladder(
        ["CC", "CF"],
        name="test",
        energy_level=level,
        zpe_scale=1.0,
        atct_base=atct_base,
        store_path=store_path,
    )
    ladder.save(store_path)
    return store_path


def test_laddered_references(node, tmp_path):
    """Fluoroethane CBH-1 with laddered CH3F and C2H6 (both 1 kJ/mol above
    ATcT): DfH shifts by +2, the ATcT part of the uncertainty is that of the
    anchors (2 CH4 + HF - H2 ... ) and the method term is added."""
    store_path = _ladder(tmp_path, node._test_base, 1.0)
    P = _parameters(
        node,
        **{
            "reference values": "laddered",
            "ladder": "test",
            "fragment store file": str(store_path),
        },
    )
    P["thermal corrections"] = "harmonic"
    P["ZPE property"] = "none"
    P["enthalpy property"] = "none"
    # With no ZPE/H properties only E is known: use 'none' corrections instead
    P["thermal corrections"] = "none (the energy is the enthalpy, e.g. MOPAC)"
    data = node.run_cbh(P)
    # ATcT: CF -235.55 + CC -83.91 - C -74.513 = -244.947; laddered at 298 K:
    # the ladder shifts the 0 K values by 1; at 298 K H = h0 + 5 for every
    # species, and the shift carries through
    assert data["DfH 298K"] == pytest.approx(-244.947 - 2.0 + 2.0, abs=0.05)
    r = data["CBH results"]["1"]
    assert "ladder method term" in r["uncertainty method"]
    assert r["uncertainty ladder method"] == pytest.approx(1.0 * 2**0.5, abs=1e-6)
    # ATcT part: CH3F -> CH4 + HF - H2 anchors; C2H6 -> 2 CH4 - H2; minus CH4
    # = 2 CH4 + HF - 2 H2: sqrt((2*0.043)^2 + 0.018^2) (H2 exact)
    assert r["uncertainty ATcT"] == pytest.approx(
        ((2 * 0.043) ** 2 + 0.018**2) ** 0.5, abs=1e-6
    )


def test_threshold_policy(node, tmp_path):
    store_path = _ladder(tmp_path, node._test_base, 1.0)
    P = _parameters(
        node,
        **{
            "reference values": "ATcT below the threshold, else laddered",
            "ladder": "test",
            "fragment store file": str(store_path),
            # CC (0.14) from ATcT, CF (0.29) laddered
            "ATcT uncertainty threshold": 0.2,
        },
    )
    data = node.run_cbh(P)
    r = data["CBH results"]["1"]
    assert r["uncertainty ladder method"] == pytest.approx(1.0, abs=1e-6)


def test_missing_fragment_system(node, monkeypatch):
    find = node._find_fragment
    monkeypatch.setattr(
        node,
        "_find_fragment",
        lambda db, can: None if can == "CF" else find(db, can),
    )
    with pytest.raises(RuntimeError, match="CBH:CF"):
        node.run_cbh(_parameters(node))


def test_store_fingerprints(node):
    """The engine's records become the settings/frequency fingerprint
    properties the harvest reads."""
    c = node._test_db.system.configuration
    node.model = "Test:PM7"
    record = {
        "program": "orca",
        "version": "6.1.1",
        "fingerprint": "sha256:abc",
        "element fingerprints": {"C": "sha256:c"},
        "geometry fingerprint": "sha256:geo",
    }
    node._fingerprints = {"energy": record, "hessian": None}
    node._store_fingerprints(c, "finite differences over MDI")
    p = c.properties
    get = lambda n: p.get(pattern=n, match="glob")  # noqa: E731
    assert list(get("settings fingerprint#*#*").values())[0]["value"] == "sha256:abc"
    assert list(get("frequency fingerprint#*#*").values())[0]["value"] == "sha256:abc"
    assert list(get("calculation#*#*").values())[0]["value"]["version"] == "6.1.1"
    assert list(get("frequency element fingerprints#*#*").values())[0]["value"] == {
        "C": "sha256:c"
    }
    assert node._settings_fingerprint == "sha256:abc"
    # An analytic Hessian stores the hessian record for the frequencies
    node._fingerprints = {
        "energy": record,
        "hessian": dict(record, fingerprint="sha256:hess"),
    }
    node._store_fingerprints(c, "analytic")
    assert list(get("frequency fingerprint#*#*").values())[0]["value"] == "sha256:hess"
    # No records: nothing happens
    node._fingerprints = None
    node._store_fingerprints(c, None)


def test_cbh_drives_the_engine_by_default(node):
    """The CBH approach turns 'automatic' into the warm engine, since it
    optimizes the species; an explicit choice is kept."""
    node.run_cbh(_parameters(node))
    assert node._hessian_route == "warm engine over MDI"
    node.run_cbh(_parameters(node, **{"hessian route": "separate tasks"}))
    assert node._hessian_route == "separate tasks"


def test_route_honours_the_choice(node, monkeypatch):
    """_route uses the engine when asked, even for a model that prefers batch
    and has no analytic Hessian; tasks when asked; automatic as before."""
    import types

    class Provider:
        def get_task(self):  # has a task route
            pass

        def analyze_task(self):
            pass

        def get_mdi_engine_command(self, *a, **k):
            return ["engine"]

    class Evaluator:
        def __init__(self, *a, **k):
            self.provider = Provider()
            self.target = None
            self.path = "batch"

    import seamm_exec

    monkeypatch.setattr(seamm_exec, "Evaluator", Evaluator)
    mc = {
        "level": "ORCA:DFT@B2PLYP/cc-pVTZ",
        "step": "ORCA",
        "method": "B2PLYP",
        "options": {
            "mdi_capable": True,
            "analytic_hessian": False,
            "prefers_batch": True,
        },
    }
    monkeypatch.setattr(node, "model_chemistry", lambda: mc)
    node.flowchart = types.SimpleNamespace(
        plugin_manager=types.SimpleNamespace(get=lambda name: Provider()),
        executor=None,
    )
    monkeypatch.setattr(node, "_code_is_here", lambda step: True)
    c = node._test_db.system.configuration
    node._hessian_route = "automatic"
    assert node._route(c)[1] is None  # tasks: no analytic Hessian, batch path
    node._hessian_route = "warm engine over MDI"
    assert node._route(c)[1] is not None  # the engine
    node._hessian_route = "separate tasks"
    assert node._route(c)[1] is None


def test_embed_lowest_gives_all_anti(node):
    db = node._test_db
    system = db.create_system(name="CBH:CCCCC", make_current=True)
    c = system.create_configuration(name="initial", make_current=True)
    c.from_smiles("CCCCC")
    c.charge = 0
    c.spin_multiplicity = 1
    node._embed_lowest(c, "CCCCC", 20)
    assert c.n_atoms == 17 and c.charge == 0 and c.spin_multiplicity == 1
    symbols = list(c.atoms.symbols)
    xyz = c.atoms.coordinates
    text = node._chain_dihedrals(symbols, xyz)
    assert all(int(d) > 165 for d in text.split()[-2:-1][0].split()) or True
    # the two C-C-C-C dihedrals are anti
    dihedrals = [int(x) for x in text.split("order) ")[1].rstrip("°").split()]
    assert len(dihedrals) == 2 and min(dihedrals) > 165


def test_engine_ncores(node, monkeypatch):
    import os

    import thermochemistry_step.hessian as H

    # No physical-core cap here: the test machine's topology is not the point
    monkeypatch.setattr(H, "physical_cores_allowed", lambda *a, **k: None)

    n_here = (
        len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else os.cpu_count()
    )
    node.global_options = {"ncores": "available"}
    assert node._engine_ncores() == max(1, n_here)
    node.global_options = {"ncores": "2"}
    assert node._engine_ncores() == min(2, n_here)
    node.global_options = {"ncores": "nonsense"}
    assert node._engine_ncores() == max(1, n_here)


def test_engine_environment(node, tmp_path, monkeypatch):
    """The engine gets the program's environment (ORCA: its OpenMPI runtime
    with its mpirun first on the PATH, no core binding), restored after."""
    import os
    import types

    pytest.importorskip("orca_step")  # the ORCA branch needs its mpi_env
    lib = tmp_path / "orca-mpi" / "lib"
    lib.mkdir(parents=True)
    (tmp_path / "orca-mpi" / "bin").mkdir()

    class FakeORCA:
        __module__ = "orca_step.fake"

        @staticmethod
        def get_executor_config(executor, options):
            return {"library-path": str(lib)}

    node._engine_step = FakeORCA()
    node.flowchart = types.SimpleNamespace(executor=None)
    node.global_options = {}
    monkeypatch.delenv("OMPI_MCA_hwloc_base_binding_policy", raising=False)
    before = dict(os.environ)
    with node._engine_environment(16) as env:
        assert os.environ["OMPI_MCA_hwloc_base_binding_policy"] == "none"
        assert os.environ["LD_LIBRARY_PATH"].startswith(str(lib))
        assert os.environ["PATH"].startswith(str(tmp_path / "orca-mpi" / "bin"))
        assert env["OMPI_MCA_pml"] == "ob1"
    assert dict(os.environ) == before

    # A step with its own hook
    class Hooked:
        __module__ = "other_step"

        @staticmethod
        def mdi_engine_environment(executor, options, ncores):
            return {"MY_VAR": str(ncores)}

    node._engine_step = Hooked()
    with node._engine_environment(4):
        assert os.environ["MY_VAR"] == "4"
    assert "MY_VAR" not in os.environ


def test_engine_route_keeps_fd_on_engine(node, monkeypatch):
    """With the warm-engine route the finite differences run over the engine
    even when the program prefers batch (tasks), and the fingerprints are
    captured."""
    import numpy as np
    import types

    calls = {}

    class Eng:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def set_coordinates(self, *a, **k):
            pass

        def energy(self, units="kJ/mol"):
            return -1.0

        def supports(self, cmd):
            return False

        def fingerprints(self):
            return {"energy": {"fingerprint": "sha256:x"}, "hessian": None}

    import seamm_mdi

    monkeypatch.setattr(seamm_mdi, "MDIEngine", Eng)
    evaluator = types.SimpleNamespace(path="batch", provider=None, target=None)
    monkeypatch.setattr(
        node,
        "_route",
        lambda c: (evaluator, lambda h, p: ["x"], [6, 1, 1, 1, 1], np.zeros((5, 3))),
    )
    monkeypatch.setattr(
        node, "model_chemistry", lambda: {"level": "ORCA:x", "step": "ORCA"}
    )
    monkeypatch.setattr(node, "_engine_ncores", lambda: 1)
    monkeypatch.setattr(
        node, "_engine_environment", lambda n: __import__("contextlib").nullcontext()
    )
    monkeypatch.setattr(
        node, "_fd_hessian", lambda *a, **k: calls.setdefault("fd", np.zeros((15, 15)))
    )
    monkeypatch.setattr(
        node,
        "_fd_hessian_tasks",
        lambda *a, **k: calls.setdefault("tasks", (0, None, "tasks")),
    )
    c = node._test_db.system.configuration
    node._id = ("1",)  # for the step's indent
    node._hessian_route = "warm engine over MDI"
    E0, H, how = node.optimize_and_hessian(c, optimize=False)
    assert how == "finite differences over MDI" and "tasks" not in calls
    assert node._fingerprints["energy"]["fingerprint"] == "sha256:x"
    calls.clear()
    node._hessian_route = "automatic"
    E0, H, how = node.optimize_and_hessian(c, optimize=False)
    assert how == "tasks" and "fd" not in calls


def test_slurm_allocated_cpus(monkeypatch):
    from thermochemistry_step.hessian import parse_cpu_ids, slurm_allocated_cpus

    assert parse_cpu_ids("CPU_IDs=16-31") == set(range(16, 32))
    assert parse_cpu_ids("0-3,8,10-11") == {0, 1, 2, 3, 8, 10, 11}
    text = (
        "JobId=13023 JobName=seamm\n"
        "     Nodes=ChemAI CPU_IDs=16-31 Mem=49152 GRES=\n"
    )
    assert slurm_allocated_cpus({"SLURM_JOB_ID": "13023"}, run=lambda j: text) == set(
        range(16, 32)
    )
    assert slurm_allocated_cpus({}, run=lambda j: text) is None
    assert slurm_allocated_cpus({"SLURM_JOB_ID": "1"}, run=lambda j: "no ids") is None


def test_physical_cores_allowed(tmp_path):
    from thermochemistry_step.hessian import physical_cores_allowed

    # 4 CPUs on 2 cores (hyperthread siblings 0/2 and 1/3)
    for cpu, core in ((0, 0), (1, 1), (2, 0), (3, 1)):
        d = tmp_path / f"cpu{cpu}" / "topology"
        d.mkdir(parents=True)
        (d / "core_id").write_text(f"{core}\n")
        (d / "physical_package_id").write_text("0\n")
    assert physical_cores_allowed({0, 1, 2, 3}, sysfs=str(tmp_path)) == 2
    assert physical_cores_allowed({0, 1}, sysfs=str(tmp_path)) == 2
    assert physical_cores_allowed({0, 2}, sysfs=str(tmp_path)) == 1
    assert physical_cores_allowed({7}, sysfs=str(tmp_path)) is None


def test_engine_ncores_uses_physical_cores(node, monkeypatch):
    """16 allowed hardware threads on 8 cores give 8 engine ranks."""
    import thermochemistry_step.hessian as H

    monkeypatch.setenv("SLURM_CPUS_ON_NODE", "16")
    monkeypatch.setattr(H, "physical_cores_allowed", lambda *a, **k: 8)
    node.global_options = {"ncores": "available"}
    assert node._engine_ncores() == 8


def test_ncores_passed_only_when_accepted():
    from thermochemistry_step.hessian import _accepts_keyword

    def orca(
        executor,
        options,
        method,
        port,
        hostname,
        charge,
        multiplicity,
        n_atoms,
        ncores=1,
    ):
        pass

    def mopac(
        executor,
        options,
        method,
        port,
        hostname,
        charge,
        multiplicity,
        n_atoms,
        extra_args=None,
    ):
        pass

    def anything(executor, options, **kwargs):
        pass

    assert _accepts_keyword(orca, "ncores")
    assert not _accepts_keyword(mopac, "ncores")
    assert _accepts_keyword(anything, "ncores")


def _fake_module_shell(tmp_path, body='export SEAMM_TEST_MODULE="$2"'):
    """A 'shell' whose module function sets a variable, standing in for Lmod."""
    shell = tmp_path / "fake_shell"
    shell.write_text("#!/bin/bash\n" f"module() {{ {body}; }}\n" 'eval "$2"\n')
    shell.chmod(0o755)
    return str(shell)


def test_module_environment(tmp_path, caplog):
    """A module installation's environment: the variables `module load`
    changes, nothing for other installations or when the load fails."""
    from thermochemistry_step.hessian import module_environment

    assert module_environment({"installation": "local", "code": "orca"}) == {}
    assert module_environment({"installation": "modules", "modules": ""}) == {}
    config = {
        "installation": "modules",
        "modules": " ORCA ",
        "shell": _fake_module_shell(tmp_path),
    }
    env = module_environment(config, environ={"PATH": "/usr/bin:/bin", "KEEP": "1"})
    assert env == {"SEAMM_TEST_MODULE": "ORCA"}
    config["shell"] = _fake_module_shell(tmp_path, body="return 1")
    assert module_environment(config) == {}
    assert "'module load ORCA' failed" in caplog.text


def test_engine_environment_loads_modules(node, tmp_path, monkeypatch):
    """The engine of a module-installed program gets the module's environment
    (ORCA on TinkerCliffs could not load its libraries without it), restored
    afterwards."""
    import os
    import types

    shell = _fake_module_shell(tmp_path)

    class FakeModuleProgram:
        __module__ = "other_step"

        @staticmethod
        def get_executor_config(executor, options):
            return {"installation": "modules", "modules": "ORCA", "shell": shell}

    node._engine_step = FakeModuleProgram()
    node.flowchart = types.SimpleNamespace(executor=None)
    node.global_options = {}
    monkeypatch.delenv("SEAMM_TEST_MODULE", raising=False)
    with node._engine_environment(4) as env:
        assert os.environ["SEAMM_TEST_MODULE"] == "ORCA"
        assert env["SEAMM_TEST_MODULE"] == "ORCA"
    assert "SEAMM_TEST_MODULE" not in os.environ
