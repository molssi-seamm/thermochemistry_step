# -*- coding: utf-8 -*-

"""The species' optimizer: ASE's BFGS (the default) or the SEAMM optimizer, over
the Model Chemistry's engine (seamm_exec's driver) or through the subflowchart."""

import numpy as np
import pytest

from seamm_util import Q_

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")


def _P(**changes):
    P = {
        "optimize fragments": "yes",
        "fragment max steps": 50,
        "fragment convergence": Q_(0.005, "eV/Å"),
        "thermal corrections": "none",
        "fragment optimizer": "ASE BFGS",
        "optimizer cheap level": "xTB:SQM@GFN2-xTB",
    }
    P.update(changes)
    return P


def test_optimizer_kwargs():
    from thermochemistry_step.cbh_approach import CBHApproach

    kwargs = CBHApproach._optimizer_kwargs
    # ASE BFGS is the driver's default: nothing to pass (and older flowcharts
    # without the parameter get it too)
    assert kwargs(_P()) == {}
    P = _P()
    del P["fragment optimizer"]
    assert kwargs(P) == {}
    # the SEAMM optimizer with a cheap level, or with Lindh's model Hessian
    assert kwargs(_P(**{"fragment optimizer": "SEAMM optimizer"})) == {
        "optimizer": "SEAMM",
        "optimizer_options": {},
        "cheap_level": "xTB:SQM@GFN2-xTB",
    }
    assert kwargs(
        _P(**{"fragment optimizer": "SEAMM optimizer", "optimizer cheap level": "none"})
    ) == {
        "optimizer": "SEAMM",
        "optimizer_options": {"hessian": "lindh"},
        "cheap_level": None,
    }


@pytest.mark.parametrize("choice", ["ASE BFGS", "SEAMM optimizer"])
def test_the_choice_reaches_the_driver(node, choice):
    """A species computed through the Model Chemistry passes the choice on."""
    from molsystem import SystemDB

    seen = []

    def optimize_and_energy(configuration, **kwargs):
        seen.append(kwargs)
        node._last_optimization = {"converged": True}
        return -100.0

    node.optimize_and_energy = optimize_and_energy
    node._store_species_energy = lambda configuration, E: None
    node._store_fingerprints = lambda *args: None
    db = SystemDB(filename="file:choice?mode=memory&cache=shared")
    configuration = db.create_system(name="CCO").create_configuration()
    configuration.from_smiles("CCO")

    result = node._compute_species_mc(
        _P(**{"fragment optimizer": choice}), "CCO", configuration
    )
    assert result["converged"] is True
    assert len(seen) == 1
    if choice == "ASE BFGS":
        assert "optimizer" not in seen[0]
    else:
        assert seen[0]["optimizer"] == "SEAMM"
        assert seen[0]["cheap_level"] == "xTB:SQM@GFN2-xTB"
    # the continuation is still the driver's
    assert seen[0]["continuations"] == 1


def test_seamm_optimizer_through_the_subflowchart(node):
    """Without an engine (the subflowchart through ASE's calculator) the SEAMM
    optimizer runs with Lindh's model Hessian and leaves the structure and its
    energy, as ASE's BFGS does."""
    from ase.build import molecule
    from ase.calculators.emt import EMT

    atoms = molecule("CH3OH")
    atoms.rattle(0.05, seed=3)
    atoms.calc = EMT()

    class Coordinates:
        def set_coordinates(self, xyz, fractionals=False):
            self.xyz = np.array(xyz)

    class Configuration:
        atoms = Coordinates()

    configuration = Configuration()
    node.ase_atoms = lambda configuration: atoms
    said = []
    node._say_continuing = lambda can, n: said.append(n)

    fmax = 0.005
    converged = node._optimize(
        configuration,
        _P(**{"fragment optimizer": "SEAMM optimizer", "fragment max steps": 200}),
        "CO",
        fmax=fmax,
    )
    assert converged
    assert said == []
    assert np.allclose(configuration.atoms.xyz, atoms.positions)
    assert np.max(np.linalg.norm(atoms.get_forces(), axis=1)) <= fmax
    assert node._last_energy == pytest.approx(
        Q_(atoms.get_potential_energy(), "eV").m_as("kJ/mol")
    )
