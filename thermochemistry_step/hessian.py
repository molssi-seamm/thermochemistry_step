# -*- coding: utf-8 -*-

"""The energy and Hessian of a structure from the flowchart's Model Chemistry.

A thin mixin over :class:`seamm_exec.hessian.HessianDriver`, which holds the
route choice (warm MDI engine or tasks; analytic or finite differences), the
engine's environment and cores, and the optimization over the engine. The
driver moved to seamm_exec on 2026-11-05 so that the Thermochemistry, Normal
Mode Sampling and Vibrational Analysis steps share one copy.
"""

import logging

import seamm_util.printing as printing

# The helpers moved with the driver; imported here for older callers.
from seamm_exec.hessian import (  # noqa: F401
    HessianDriver,
    _accepts_keyword,
    module_environment,
    modules_to_load,
    parse_cpu_ids,
    physical_cores_allowed,
    slurm_allocated_cpus,
)

logger = logging.getLogger(__name__)
printer = printing.getPrinter("Thermochemistry")


class ModelChemistryHessian:
    """Mixin: the energy and Cartesian Hessian from the Model Chemistry.

    The route is the step's ``_hessian_route`` ("automatic", "warm engine over
    MDI" or "separate tasks"); the engine's settings fingerprints from the last
    calculation are in ``_fingerprints``.
    """

    def model_chemistry(self):
        """The ``_model_chemistry`` wrapper: a preset's own (see the CBH
        approach's presets), else the one published upstream, or None."""
        preset = getattr(self, "_preset_mc", None)
        if preset is not None:
            return preset
        if not self.variable_exists("_model_chemistry"):
            return None
        return self.get_variable("_model_chemistry")

    def _hessian_driver(self):
        return HessianDriver(
            self,
            self.model_chemistry(),
            route=getattr(self, "_hessian_route", "automatic") or "automatic",
            who="Thermochemistry",
            hint="or put the calculation in this step's subflowchart",
            print_to=printer,
        )

    def model_chemistry_energy(self, configuration):
        """The energy (kJ/mol) of the configuration with the Model Chemistry."""
        return self._hessian_driver().energy(configuration)

    def energy_and_hessian(self, configuration, delta_bohr=0.01):
        """The energy (kJ/mol), the Cartesian Hessian (hartree/bohr², (3N, 3N))
        and how it was obtained: "analytic", "finite differences over MDI" or
        "finite differences as tasks"."""
        return self.optimize_and_hessian(configuration, delta_bohr, optimize=False)

    def optimize_and_energy(
        self,
        configuration,
        optimize=True,
        fmax=0.005,
        max_steps=200,
        continuations=0,
        on_continue=None,
        optimizer="ASE BFGS",
        optimizer_options=None,
        cheap_level=None,
    ):
        """Optionally optimize the structure over the warm engine, then its
        energy (kJ/mol). With no engine (tasks) the energy only. With
        `continuations` an unconverged optimization continues in the same engine
        session, calling ``on_continue(max_steps)`` first. `optimizer` ("ASE BFGS"
        or "SEAMM"), `optimizer_options` and `cheap_level` are as in seamm_exec's
        HessianDriver.optimize."""
        driver = self._hessian_driver()
        try:
            return driver.optimize_and_energy(
                configuration,
                optimize,
                fmax,
                max_steps,
                continuations=continuations,
                on_continue=on_continue,
                optimizer=optimizer,
                optimizer_options=optimizer_options,
                cheap_level=cheap_level,
            )
        finally:
            self._fingerprints = driver.fingerprints
            self._last_optimization = getattr(driver, "last_optimization", None)
            # (<S**2>, S(S+1)) at the end of the optimization, when the engine
            # gives them (ORCA's), for the species' properties
            self._spin = getattr(driver, "spin", None)

    def optimize_and_hessian(
        self,
        configuration,
        delta_bohr=0.01,
        optimize=True,
        fmax=0.005,
        max_steps=200,
        continuations=0,
        on_continue=None,
        optimizer="ASE BFGS",
        optimizer_options=None,
        cheap_level=None,
        finite_differences=False,
        check=False,
    ):
        """Optionally optimize the structure (ASE BFGS or the SEAMM optimizer over
        the warm engine), then the energy and Hessian, in one engine session.
        With `finite_differences` the Hessian is central differences of the
        engine's forces even when it has an analytic one. With
        `continuations` an unconverged optimization continues in the same
        session, so the Hessian is computed once, at the final structure.
        With `check` a finite-difference Hessian is also tested against
        differences of the gradient (seamm_exec raises ``HessianRefused`` if it,
        or any of its points, fails). ``self._fd_info`` and
        ``self._hessian_check`` say how it was made and checked.

        Returns
        -------
        (float, numpy.ndarray, str)
            As `energy_and_hessian`.
        """
        driver = self._hessian_driver()
        try:
            return driver.optimize_and_hessian(
                configuration,
                delta_bohr,
                optimize,
                fmax,
                max_steps,
                continuations=continuations,
                on_continue=on_continue,
                optimizer=optimizer,
                optimizer_options=optimizer_options,
                cheap_level=cheap_level,
                finite_differences=finite_differences,
                check=check,
            )
        finally:
            self._fd_info = getattr(driver, "fd_info", None)
            self._hessian_check = getattr(driver, "last_check", None)
            self._fingerprints = driver.fingerprints
            self._last_optimization = getattr(driver, "last_optimization", None)
            # (<S**2>, S(S+1)) at the end of the optimization, when the engine
            # gives them (ORCA's), for the species' properties
            self._spin = getattr(driver, "spin", None)
