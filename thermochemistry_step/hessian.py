# -*- coding: utf-8 -*-

"""The energy and Hessian of a structure from the flowchart's Model Chemistry.

The Model Chemistry step upstream publishes ``_model_chemistry``; the program
behind it is driven through ``seamm_exec.Evaluator``, either as a warm MDI
engine (``seamm_mdi.MDIEngine``) or as separate calculations (tasks) where the
job runs them. The route follows the Normal Mode Sampling step (design of
2026-10-03):

* on a queue target the finite difference runs as tasks on the cluster, with
  no local engine started even to ask;
* otherwise the analytic Hessian over MDI when the program offers one (the
  engine advertises the custom ``<HESSIAN`` command only when the method has
  a genuine analytic Hessian), else central finite differences of the forces
  over the warm engine, or as tasks when the evaluator chose that path;
* with no MDI engine, the finite difference as tasks.

This mirrors ``normal_mode_sampling_step``'s driver, with the energy of the
undisplaced structure collected as well. It belongs in a shared library once a
third consumer appears.
"""

import logging

import numpy as np

from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
printer = printing.getPrinter("Thermochemistry")


class _MDICalculator:
    """A minimal ASE calculator over a warm MDI engine, for the optimization."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, engine):
        self.engine = engine
        self.results = {}
        self.n_calls = 0
        self._positions = None

    def _update(self, atoms):
        if self._positions is None or not np.array_equal(
            self._positions, atoms.positions
        ):
            self.engine.set_coordinates(atoms.positions, units="Å")
            self.results = {
                "energy": self.engine.energy(units="eV"),
                "forces": np.asarray(self.engine.forces(units="eV/Å"), dtype=float),
            }
            self._positions = np.array(atoms.positions)
            self.n_calls += 1

    def get_potential_energy(self, atoms=None, force_consistent=False):
        self._update(atoms)
        return self.results["energy"]

    def get_forces(self, atoms=None):
        self._update(atoms)
        return self.results["forces"]

    def calculation_required(self, atoms, properties):
        return self._positions is None or not np.array_equal(
            self._positions, atoms.positions
        )


class ModelChemistryHessian:
    """Mixin: the energy and Cartesian Hessian from the Model Chemistry."""

    def model_chemistry(self):
        """The ``_model_chemistry`` wrapper published upstream, or None."""
        if not self.variable_exists("_model_chemistry"):
            return None
        return self.get_variable("_model_chemistry")

    def _route(self, configuration):
        """How to drive the Model Chemistry for this configuration.

        Returns
        -------
        (evaluator, build_argv or None, elements, coords_ang)
            ``build_argv`` is None when the finite difference must run as
            tasks (no MDI engine, a queue target, or the program not here).
        """
        mc = self.model_chemistry()
        if mc is None:
            raise ValueError(
                "Thermochemistry needs a Model Chemistry: add a 'Model Chemistry' "
                "step before it to choose the method, or put the calculation in "
                "this step's subflowchart."
            )
        if configuration.periodicity != 0:
            raise ValueError(
                "The thermochemistry of a periodic system is not supported."
            )

        options = mc["options"]
        from seamm_exec import Evaluator

        try:
            evaluator = Evaluator(self, mc, properties=("energy", "gradients"))
        except ValueError:
            raise ValueError(
                f"The model chemistry '{mc['level']}' can be evaluated neither over "
                "MDI nor as separate calculations, so it cannot give the Hessian. "
                "Choose a model chemistry with an MDI engine or with tasks (e.g. "
                "MOPAC, xTB or ORCA)."
            )

        step = self.flowchart.plugin_manager.get(mc["step"])
        method = options.get("mdi_method_arg") or mc["method"]
        basis = options.get("mdi_basis_arg")

        def build_argv(hostname, port):
            kwargs = {}
            if basis is not None:
                kwargs["basis"] = basis
            return step.get_mdi_engine_command(
                self.flowchart.executor,
                self.global_options,
                method=method,
                port=port,
                hostname=hostname,
                charge=configuration.charge,
                multiplicity=configuration.spin_multiplicity,
                n_atoms=configuration.n_atoms,
                **kwargs,
            )

        elements = list(configuration.atoms.atomic_numbers)
        coords_ang = np.asarray(
            configuration.atoms.get_coordinates(fractionals=False, as_array=True),
            dtype=float,
        )

        has_tasks = hasattr(evaluator.provider, "get_task")
        target = evaluator.target
        on_queue = target is not None and getattr(target, "tasks", None) in (
            "queue",
            "taskserver",
        )
        use_tasks = False
        if on_queue and has_tasks:
            use_tasks = True
        elif not options.get("mdi_capable", False):
            use_tasks = True
        elif options.get("analytic_hessian") is False and evaluator.path == "batch":
            use_tasks = True
        elif has_tasks and not self._code_is_here(step):
            printer.important(
                __(
                    f"{mc['level']} is not installed on this machine, so the "
                    "Hessian is computed by finite differences, as separate "
                    "calculations where the job runs its tasks.",
                    indent=self.indent + 4 * " ",
                )
            )
            use_tasks = True
        return evaluator, (None if use_tasks else build_argv), elements, coords_ang

    def model_chemistry_energy(self, configuration):
        """The energy (kJ/mol) of the configuration with the Model Chemistry."""
        from seamm_exec import Evaluator

        mc = self.model_chemistry()
        if mc is None:
            raise ValueError("Thermochemistry needs a Model Chemistry step before it.")
        with Evaluator(self, mc, properties=("energy",)) as evaluator:
            evaluator.submit(configuration, key="energy")
            for result in evaluator.results():
                if not result.ok:
                    raise RuntimeError(f"The energy failed: {result.reason}")
                return float(result.energy)
        raise RuntimeError("No energy was returned by the Model Chemistry")

    def energy_and_hessian(self, configuration, delta_bohr=0.01):
        """The energy (kJ/mol) and the Cartesian Hessian (hartree/bohr²,
        (3N, 3N)) of the configuration with the Model Chemistry.

        Returns
        -------
        (float, numpy.ndarray, str)
            The energy, the Hessian, and how the Hessian was obtained:
            "analytic", "finite differences over MDI" or "finite differences
            as tasks".
        """
        return self.optimize_and_hessian(configuration, delta_bohr, optimize=False)

    def optimize_and_energy(
        self, configuration, optimize=True, fmax=0.005, max_steps=200
    ):
        """Optionally optimize the structure over the warm engine, then its
        energy (kJ/mol). With no engine (tasks) the energy only."""
        evaluator, build_argv, elements, coords_ang = self._route(configuration)
        if build_argv is None or not optimize or len(elements) == 1:
            return self.model_chemistry_energy(configuration)

        from seamm_mdi import MDIEngine

        with MDIEngine(build_argv, elements, name="SEAMM", logger=self.logger) as eng:
            coords_ang = self._optimize_over_engine(
                eng, configuration, elements, coords_ang, fmax, max_steps
            )
            eng.set_coordinates(coords_ang, units="Å")
            return eng.energy(units="kJ/mol")

    def _optimize_over_engine(
        self, eng, configuration, elements, coords_ang, fmax, max_steps
    ):
        """ASE BFGS over the warm engine; updates the configuration."""
        from ase import Atoms
        from ase.optimize import BFGS

        atoms = Atoms(numbers=elements, positions=coords_ang)
        atoms.calc = _MDICalculator(eng)
        converged = BFGS(atoms, logfile=None).run(fmax=fmax, steps=max_steps)
        coords_ang = np.asarray(atoms.positions, dtype=float)
        configuration.atoms.set_coordinates(coords_ang, fractionals=False)
        printer.important(
            __(
                f"The structure was optimized in {atoms.calc.n_calls} energy "
                "evaluations"
                + ("." if converged else f" without converging in {max_steps} steps."),
                indent=self.indent + 4 * " ",
            )
        )
        return coords_ang

    def optimize_and_hessian(
        self,
        configuration,
        delta_bohr=0.01,
        optimize=True,
        fmax=0.005,
        max_steps=200,
    ):
        """Optionally optimize the structure (ASE BFGS over the warm engine),
        then the energy and Hessian, in one engine session.

        Parameters
        ----------
        configuration : molsystem._Configuration
            Its coordinates are updated by the optimization.
        delta_bohr : float
            The finite-difference step.
        optimize : bool
        fmax : float
            The convergence criterion for the optimization, eV/Å.
        max_steps : int

        Returns
        -------
        (float, numpy.ndarray, str)
            As `energy_and_hessian`.
        """
        evaluator, build_argv, elements, coords_ang = self._route(configuration)
        n = len(elements)
        if build_argv is None:
            if optimize:
                printer.important(
                    __(
                        "The structure is not optimized: the Model Chemistry runs as "
                        "separate tasks, with no engine to drive the optimization.",
                        indent=self.indent + 4 * " ",
                    )
                )
            return self._fd_hessian_tasks(
                evaluator, configuration, coords_ang, delta_bohr
            )

        from seamm_mdi import MDIEngine  # only for the MDI routes

        with MDIEngine(build_argv, elements, name="SEAMM", logger=self.logger) as eng:
            if optimize and n > 1:
                coords_ang = self._optimize_over_engine(
                    eng, configuration, elements, coords_ang, fmax, max_steps
                )
            eng.set_coordinates(coords_ang, units="Å")
            E0 = eng.energy(units="kJ/mol")
            if eng.supports("<HESSIAN"):
                H = np.asarray(eng.hessian(), dtype=float).reshape(3 * n, 3 * n)
                return E0, H, "analytic"
            if evaluator.path != "batch":
                H = self._fd_hessian(eng, coords_ang, n, delta_bohr)
                return E0, H, "finite differences over MDI"
        return self._fd_hessian_tasks(evaluator, configuration, coords_ang, delta_bohr)

    def _code_is_here(self, step):
        """Whether the program is installed where this step runs."""
        check = getattr(step, "get_executor_config", None)
        if check is None:
            return True
        try:
            check(self.flowchart.executor, self.global_options)
        except Exception:
            return False
        return True

    def _fd_hessian_tasks(self, evaluator, configuration, coords_ang, delta_bohr=0.01):
        """The energy and the central-difference Hessian (hartree/bohr²) from
        the gradients of the 6N displaced structures, computed together as tasks
        (concurrently, or on the job's queue; a rerun keeps the finished ones)."""
        from seamm_exec import Geometry

        elements = list(configuration.atoms.atomic_numbers)
        x0 = np.asarray(coords_ang, dtype=float).reshape(-1)
        delta = Q_(delta_bohr, "bohr").m_as("Å")
        ndof = x0.size
        prefix = f"c{getattr(configuration, 'id', 0)}-fd"
        evaluator.submit(
            Geometry(
                elements,
                x0.reshape(-1, 3),
                configuration.charge,
                configuration.spin_multiplicity,
            ),
            key=f"{prefix}eq",
        )
        for j in range(ndof):
            for sign, tag in ((1.0, "p"), (-1.0, "m")):
                x = x0.copy()
                x[j] += sign * delta
                evaluator.submit(
                    Geometry(
                        elements,
                        x.reshape(-1, 3),
                        configuration.charge,
                        configuration.spin_multiplicity,
                    ),
                    key=f"{prefix}{j:04d}{tag}",
                )
        gradients = {}
        E0 = None
        for result in evaluator.results():
            if not result.ok:
                raise RuntimeError(
                    f"A structure for the Hessian failed ({result.key}): "
                    f"{result.reason}"
                )
            if result.key == f"{prefix}eq":
                E0 = float(result.energy)
                continue
            gradients[result.key] = (
                Q_(np.asarray(result.gradients, dtype=float), "kJ/mol/Å")
                .m_as("hartree/bohr")
                .reshape(-1)
            )
        H = np.zeros((ndof, ndof))
        for j in range(ndof):
            gp = gradients[f"{prefix}{j:04d}p"]
            gm = gradients[f"{prefix}{j:04d}m"]
            H[:, j] = (gp - gm) / (2.0 * delta_bohr)
        return E0, 0.5 * (H + H.T), "finite differences as tasks"

    def _fd_hessian(self, engine, coords_ang, n_atoms, delta_bohr=0.01):
        """Central-difference Hessian (hartree/bohr²) from MDI forces."""
        x0 = Q_(coords_ang, "Å").m_as("bohr").reshape(-1)
        ndof = 3 * n_atoms
        H = np.zeros((ndof, ndof))
        for j in range(ndof):
            xp = x0.copy()
            xp[j] += delta_bohr
            engine.set_coordinates(xp, units="bohr")
            fp = np.asarray(engine.forces(units="hartree/bohr"), dtype=float).reshape(
                -1
            )
            xm = x0.copy()
            xm[j] -= delta_bohr
            engine.set_coordinates(xm, units="bohr")
            fm = np.asarray(engine.forces(units="hartree/bohr"), dtype=float).reshape(
                -1
            )
            # Force = -dE/dx, so d²E/dx_i dx_j = -(f_i(+) - f_i(-)) / 2 delta.
            H[:, j] = -(fp - fm) / (2.0 * delta_bohr)
        # Leave the engine at the undisplaced structure
        engine.set_coordinates(x0, units="bohr")
        return 0.5 * (H + H.T)
