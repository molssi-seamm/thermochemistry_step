# -*- coding: utf-8 -*-

"""The isodesmic (CBH) approach with ATcT references.

The enthalpy of formation of the target molecule is obtained from an
error-cancelling reaction of the connectivity-based hierarchy (CBH), whose
other species have reference enthalpies of formation in the Active
Thermochemical Tables (ATcT):

    ΔfH(target, T) = Σ_j ν_j ΔfH_j(ATcT, T) - Δ_rH(T)

with Δ_rH(T) = Σ_i ν_i H_i(T) computed with the subflowchart's engine for
every species (the target has ν = -1). The reaction, the ATcT lookups, the
covariance-propagated 95 % uncertainty and the combination are done here; the
science lives in ``seamm_thermochemistry`` (``atct`` and ``cbh``).

The energies H_i(T) come from one of three places:

* "finite-difference harmonic": the subflowchart is run on each species, first
  to optimize it (ASE BFGS) and then for the finite-difference Hessian; the
  rigid-rotor harmonic-oscillator thermochemistry of the Harmonic approach
  gives H(0 K) = E + ZPE and H(298.15 K);
* "finite-difference harmonic" with "thermal corrections" = none: the
  subflowchart's energy already is the enthalpy of formation at 298.15 K, as
  for MOPAC's NDDO methods, so only Δ_rH(298.15 K) is formed;
* "existing results (run nothing)": the species' energies are read from
  properties already stored on the systems named ``CBH:<SMILES>`` (and the
  target), e.g. by an earlier parallel Loop, so the expensive part can be
  run in parallel and at two levels of theory.

Results are reported at 0 K and 298.15 K and paired with the ATcT values at
the same temperature, never converted between them. ATcT's plain gas-phase
entries are the lowest conformer at 0 K and a conformer average at 298.15 K,
so the 0 K comparison is the cleaner one for a single-conformer calculation.
"""

from collections import Counter
import json
import logging
import math
import os
from pathlib import Path

import numpy as np
from tabulate import tabulate

from seamm_util import Q_, parse_list
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
printer = printing.getPrinter("Thermochemistry")

T_REF = 298.15

#: The presets of the CBH approach (thermo's defaults, adopted 2026-11-08): the
#: composite recipe (composite_step), and its geometry/frequency level and energy
#: level (None: the geometry level). With composite_step installed the levels come
#: from the recipe; these spellings are the fallback without it, and a test keeps
#: them equal to the recipes'.
PRESETS = {
    "quick": {
        "recipe": "SEAMM:CMP@r2SCAN-3c-HO",
        "geometry": "ORCA:CMP@R2SCAN-3C",
        "energy": None,
    },
    "default": {
        "recipe": "SEAMM:CMP@revDSD-r2SCAN-3c-HO",
        "geometry": "ORCA:CMP@R2SCAN-3C",
        "energy": "ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP",
    },
    # The energy level is chosen per target (high_level.py): canonical or DLPNO
    # CCSD(T)-F12, the same for every species of the reaction
    "high": {
        "recipe": None,
        "geometry": "ORCA:CMP@R2SCAN-3C",
        "energy": None,  # replaced by HIGH_ENERGY below
    },
}

#: The 'high' preset's energy level, chosen when the target is known
HIGH_ENERGY = (
    "CCSD(T)-F12/cc-pVTZ-F12 (canonical, or DLPNO for a large saturated "
    "closed-shell target)"
)
PRESETS["high"]["energy"] = HIGH_ENERGY
#: The default largest canonical CCSD(T)-F12 of the 'high' preset (96 GB nodes)
CANONICAL_MAX_BF_DEFAULT = 546

#: The owner and type of SEAMM's composite recipes as model chemistries.
COMPOSITE_PREFIX = "SEAMM:CMP@"


def is_composite(level):
    """Whether a model chemistry is one of SEAMM's composite recipes."""
    return str(level or "").startswith(COMPOSITE_PREFIX)


def recipe_levels(level, data_path=None):
    """The levels of a composite recipe ``SEAMM:CMP@<recipe>`` (composite_step).

    Returns
    -------
    dict
        ``recipe`` (the level spec), ``name``, ``version``, ``fingerprint``,
        ``geometry`` (the geometry/frequency level), ``energy`` (the energy level:
        None when E0 is the geometry level's energy, the component's level when it
        is one other energy, or the composite itself, ``recipe``, when it combines
        several), ``components`` ({name: level}), ``zpe`` ((factor, source)) and
        ``object`` (the composite_step Recipe). The step runs ordinary levels, so
        the species' cached properties, the fragment store's rows and the ZPE table
        are those of the component levels; the recipe's name never keys anything.
    """
    try:
        from composite_step.evaluate import component_levels, single_source
        from composite_step.recipes import RecipeLibrary
    except ImportError:
        raise RuntimeError(
            f"The composite method {level} needs the Composite step "
            "(composite-step) installed."
        )
    name = str(level)[len(COMPOSITE_PREFIX) :]
    extra = [Path(p) / "recipes" for p in (data_path or ())]
    library = RecipeLibrary(extra_directories=extra)
    try:
        recipe = library.get(name)
    except KeyError as e:
        raise RuntimeError(f"No composite recipe for {level}: {e}")
    source = single_source(recipe)
    levels = component_levels(recipe)
    if source is None:
        # Several calculations: the composite itself is the energy level
        energy = recipe.model_chemistry
    elif source == "geometry":
        energy = None
    else:
        energy = levels[source]
    return {
        "recipe": recipe.model_chemistry,
        "name": recipe.name,
        "version": recipe.version,
        "fingerprint": recipe.fingerprint(),
        "geometry": recipe.geometry["level"],
        "energy": energy,
        "components": levels,
        "zpe": recipe.zpe_scale,
        "object": recipe,
    }


_PROTECT = {
    "default": "default",
    "none": "none",
    "aromatic rings": ("aromatic",),
    "aromatic + small rings": ("aromatic", "small rings"),
    "aromatic, small rings and charged groups": "default",
    "all rings": ("all rings", "groups"),
}


class CBHApproach:
    """Mixin for Thermochemistry: the CBH approach."""

    def describe_cbh(self, P):
        """The description of the CBH approach for the step's summary."""
        if species_only(P):
            preset = P.get("preset", "custom")
            if preset in PRESETS:
                levels = PRESETS[preset]
                how = (
                    f"with the '{preset}' preset: the geometry and frequencies with "
                    f"{levels['geometry']}, the energy with "
                    f"{levels['energy'] or levels['geometry']}"
                )
            else:
                how = "with the flowchart's Model Chemistry"
            return (
                "Compute only the system itself as a CBH species, for the fragment "
                f"store, {how}: optimized, its harmonic frequencies and its energies "
                "stored with the settings fingerprints the store's harvest reads. "
                "No CBH reactions are run and ATcT is not used."
            )
        text = (
            "Calculate the enthalpy of formation of the system with the isodesmic "
            "reactions of the connectivity-based hierarchy (CBH), using the Active "
            "Thermochemical Tables (ATcT) for the reference species."
        )
        rungs = P["cbh rungs"]
        if rungs == "all usable":
            text += (
                " Every rung from CBH-0 to CBH-4 whose reference species are all in "
                "ATcT will be used and the results compared."
            )
        else:
            text += f" The rungs {rungs} will be used where possible."
        protect = P["protect"]
        if protect == "none":
            text += " No groups will be kept whole (plain Kekulé CBH)."
        else:
            text += f" Groups kept whole: {protect}."
        mode = P["cbh energies"]
        if mode == "existing results (run nothing)":
            text += (
                " The energies of the species will be read from the properties "
                "already stored on the systems (nothing will be run)."
            )
        else:
            if P["thermal corrections"] == "harmonic":
                text += (
                    " The subflowchart will be run on the target and on each "
                    "fragment, optimizing the structure if requested, and the "
                    "finite-difference harmonic thermochemistry will give H(0 K) "
                    "and H(298.15 K)."
                )
            else:
                text += (
                    " The subflowchart's energy will be taken as the enthalpy at "
                    "298.15 K (as for MOPAC), with no thermal corrections."
                )
        if _yes(P["compare with ATcT"]):
            text += " If the target is itself in ATcT the result will be compared."
        return text

    # -----------------------------------------------------------------------

    def run_cbh(self, _P):
        """Run the CBH approach. Returns the data dictionary stored as results."""
        from seamm_thermochemistry import (
            ATcT,
            AmbiguousEntry,
            MissingEntry,
            cbh,
            reaction_uncertainty,
        )
        from seamm_thermochemistry.atct import canonical_smiles

        from .resonance import lookup_target

        _P = self.apply_preset(_P)
        self._reset_fingerprints()
        preset = _P.get("preset", "custom")
        system_db = self.get_variable("_system_db")
        target_system, target_configuration = self.get_system_configuration()
        directory = self._species_directory()
        if directory is not None and (directory / "species.sdf").exists():
            (directory / "species.sdf").unlink()
        n_atoms = target_configuration.n_atoms
        if target_configuration.bonds.n_bonds == 0 and n_atoms > 1:
            raise RuntimeError(
                "The CBH approach needs the bonds of the molecule: build it from "
                "SMILES, or read a file with bonds (SDF, MOL)."
            )
        charge = target_configuration.charge
        multiplicity = target_configuration.spin_multiplicity or 1
        smiles = self._target_smiles(target_configuration, charge, multiplicity)
        target = canonical_smiles(smiles)
        self._high_choice = None
        self._high_max_bf = _P.get(
            "maximum basis functions for canonical F12", CANONICAL_MAX_BF_DEFAULT
        )
        self._high_target = (
            smiles,
            multiplicity,
            list(target_configuration.atoms.atomic_numbers),
        )

        if species_only(_P):
            return self._run_species_only(
                _P, smiles, target, target_system, target_configuration
            )

        text = ""
        if preset in PRESETS:
            levels = PRESETS[preset]
            text += (
                f"The '{preset}' preset ({levels['recipe']}): geometries and "
                "frequencies with "
                f"{levels['geometry']}, energies with "
                f"{levels['energy'] or levels['geometry']}, the highest CBH rung "
                f"available up to CBH-{max(self._preset_rungs(_P))}.\n"
            )
        text += f"The target is {smiles}"
        if charge != 0 or multiplicity != 1:
            text += f" with charge {charge} and multiplicity {multiplicity}"
        text += "."
        # The structure's multiplicity is used; say so where the SMILES shows more
        # unpaired electrons, and where the ground state differs (a singlet carbene
        # built as a triplet, or O2 as a singlet)
        _, _, radicals, _ = cbh._composition(smiles)
        ground = cbh_multiplicity(smiles)
        if radicals > multiplicity - 1:
            text += (
                f" Its SMILES shows {radicals} unpaired electrons, but the "
                f"structure's spin multiplicity, {multiplicity}, is used"
            )
            if ground == multiplicity:
                text += f", which is the ground state of {target}."
            else:
                text += (
                    f"; check it: the ground state of {target} is taken to be "
                    f"{ground}."
                )
        elif ground != multiplicity:
            text += (
                f" Its spin multiplicity, {multiplicity}, is the structure's; the "
                f"ground state of {target} is {ground}."
            )
        text += "\n"

        # The rungs
        rungs = _P["cbh rungs"]
        if rungs in ("all usable", "default", ""):
            rungs = [0, 1, 2, 3, 4]
        else:
            rungs = sorted({int(x) for x in parse_list(str(rungs))})
        protect = _PROTECT.get(_P["protect"], _P["protect"])

        reactions = {}
        unavailable = {}
        for rung in rungs:
            try:
                reactions[rung] = cbh.cbh_reactions(
                    smiles, rungs=[rung], protect=protect
                )[rung]
            except cbh.RungUnavailable as e:
                unavailable[rung] = str(e)
        for rung, reaction in reactions.items():
            if not reaction.trivial and not reaction.balanced:
                raise RuntimeError(
                    f"The CBH-{rung} reaction does not balance: "
                    + "; ".join(reaction.balance_errors)
                    + f"\n{reaction}"
                )
        protected = next(iter(reactions.values())).protected if reactions else []
        if protected:
            text += "Groups kept whole: " + "; ".join(protected) + ".\n"

        # ATcT
        base = _P["ATcT directory"]
        base = None if base in ("default", "") else os.path.expanduser(base)
        version = _P["ATcT version"]
        version = None if version in ("latest", "") else version
        downloaded = self._ensure_atct(base, version)
        tables = ATcT(version=version, base=base)
        if downloaded:
            text += (
                f"No ATcT data were installed, so ver. {tables.version} was "
                "downloaded from atct.anl.gov.\n"
            )
        text += (
            f"Reference values from ATcT ver. {tables.version} of the "
            f"Thermochemical Network ({tables.directory}).\n"
        )

        def lookup(can):
            _, q, _, _ = cbh._composition(can)
            return tables.lookup(can, charge=q, multiplicity=cbh_multiplicity(can))

        # The reference policy: ATcT, a ladder, or ATcT below a threshold
        policy = _P.get("reference values", "ATcT")
        ladder = None
        threshold = None
        if policy != "ATcT":
            from seamm_thermochemistry import Ladder

            name = str(_P.get("ladder", "")).strip()
            if not name:
                raise RuntimeError(
                    "The reference values are to be laddered but no ladder is "
                    "named; build one with seamm-thermochemistry-ladder."
                )
            store_file = _P.get("fragment store file", "default")
            store_file = (
                None
                if store_file in ("default", "")
                else os.path.expanduser(store_file)
            )
            ladder = Ladder.load(name, store_file)
            if ladder.atct_version != tables.version:
                raise RuntimeError(
                    f"The ladder '{name}' is on ATcT ver. {ladder.atct_version} but "
                    f"the step uses ver. {tables.version}: one version per reaction."
                )
            text += (
                f"Reference values: {policy}, ladder '{name}' (energies "
                f"{ladder.data['energy_level']}, frequencies "
                f"{ladder.data['frequency_level']}, ZPE x {ladder.data['zpe_scale']})"
            )
            mu = ladder.method_uncertainty
            if mu.get("n"):
                text += (
                    f"; its method term, from {mu['n']} ladder species in ATcT, is "
                    f"{mu['rms_0K']:.2f} kJ/mol RMS at 0 K (an upper bound)"
                )
            text += ".\n"
            if policy.startswith("ATcT below"):
                threshold = _P["ATcT uncertainty threshold"]
                if hasattr(threshold, "m_as"):
                    threshold = threshold.m_as("kJ/mol")
                threshold = float(threshold)

        # The supplementary reference library, for species ATcT lacks
        # (only with ATcT as the reference values, not with a ladder)
        resolver = None
        library = {}
        supplementary = _P.get("supplementary references", "none")
        use_user = _yes(_P.get("use my own references", "no"))
        if supplementary == "none":
            use_user = False
        if policy == "ATcT" and (supplementary != "none" or use_user):
            from seamm_thermochemistry.atct import atct_directory
            from seamm_thermochemistry.references import LayeredReferences
            from seamm_thermochemistry.resolver import ReferenceResolver

            # The libraries sit beside the ATcT data the step uses
            where = atct_directory(base).parent
            resolver = ReferenceResolver(
                tables,
                supplementary=supplementary,
                use_user=use_user,
                multiplicity=cbh_multiplicity,
                layers=LayeredReferences(
                    where / "references.db", where / "references-official.db"
                ),
            )
            if not resolver.layers.layers:
                text += (
                    "No supplementary reference library is installed, so only "
                    "ATcT is used.\n"
                )
                resolver = None
            else:
                wanted = sorted(
                    {
                        can
                        for reaction in reactions.values()
                        if not reaction.trivial
                        for can in reaction.references
                    }
                )
                library, _, _, library_notes = resolver.resolve(wanted)
                library = {
                    k: v for k, v in library.items() if v.get("source") != "ATcT"
                }
                versions = resolver.layers.versions
                text += (
                    "Species ATcT lacks are taken from the reference library ("
                    + ", ".join(
                        f"{layer} {v or 'unversioned'}"
                        for layer, (v, _) in versions.items()
                    )
                    + f"; {supplementary}"
                    + ("; with your own references" if use_user else "")
                    + ").\n"
                )
                for note in library_notes:
                    text += note + "\n"

        def reference(can):
            """The reference entry under the policy: an ATcT entry (with
            atct_coefficients {id: 1}), the ladder's entry, or the reference
            library's (with no ATcT coefficients)."""
            atct_entry = None
            atct_error = None
            try:
                atct_entry = lookup(can)
            except (MissingEntry, AmbiguousEntry) as e:
                atct_error = e
            if atct_entry is not None:
                atct_entry = dict(atct_entry)
                atct_entry["atct_coefficients"] = {atct_entry["id"]: 1.0}
                atct_entry["source"] = "ATcT"
            if ladder is None:
                if atct_entry is None:
                    if can in library:
                        return library[can]
                    raise atct_error
                return atct_entry
            laddered = ladder.entry(can)
            if laddered is not None and laddered["source"] == "ATcT":
                laddered = None  # a CBH-0 anchor: the ATcT value itself
            if policy == "laddered":
                if laddered is not None:
                    return laddered
                if atct_entry is not None and _own_rung_is_zero(can):
                    return atct_entry
                raise MissingEntry(
                    f"'{can}' is not in the ladder '{ladder.name}'"
                    + (
                        f" (ATcT has it: {atct_entry['id']})"
                        if atct_entry is not None
                        else ""
                    )
                )
            # ATcT below the threshold, else laddered
            if atct_entry is not None and (
                atct_entry["exact"]
                or (atct_entry["uncertainty_298K"] or 0.0) <= threshold
            ):
                return atct_entry
            if laddered is not None:
                return laddered
            if atct_entry is not None:
                return atct_entry
            raise atct_error

        def _own_rung_is_zero(can):
            try:
                return cbh.cbh_reactions(can, rungs=[0], protect=protect)[0].trivial
            except cbh.CBHError:
                return True  # H2

        entries = {}
        usable = []
        missing = {}
        ambiguous = {}
        for rung, reaction in sorted(reactions.items()):
            if reaction.trivial:
                continue
            absent = []
            for can in reaction.references:
                if can in entries:
                    continue
                try:
                    entries[can] = reference(can)
                except MissingEntry as e:
                    if ladder is not None:
                        ambiguous[can] = str(e)
                    absent.append(can)
                except AmbiguousEntry as e:
                    ambiguous[can] = str(e)
                    absent.append(can)
            if absent:
                missing[rung] = absent
            else:
                usable.append(rung)

        target_entry = None
        target_via = None
        if _yes(_P["compare with ATcT"]):
            target_entry, target_via = lookup_target(
                tables, target, charge, multiplicity
            )

        # Report the reactions
        lines = []
        for rung in rungs:
            if rung in unavailable:
                lines.append(f"CBH-{rung}: not available. {unavailable[rung]}")
            elif reactions[rung].trivial:
                lines.append(f"CBH-{rung}: trivial, the rung reproduces the target.")
            elif rung in missing:
                lines.append(
                    f"{reactions[rung]}\n    not usable: no "
                    + ("reference value" if resolver is not None else "ATcT entry")
                    + " for "
                    + ", ".join(missing[rung])
                )
            else:
                lines.append(str(reactions[rung]))
        text += "\nThe reactions (the target has coefficient -1):\n    "
        text += "\n    ".join(lines) + "\n"
        for can, message in ambiguous.items():
            text += f"\n{message}\n"
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        if not usable:
            # A target that is itself a reference species (methane, benzene, the
            # methyl radical, CO): its rungs are trivial or need other species
            # that ATcT lacks, but its own value is the answer
            if target_entry is None:
                target_entry, target_via = lookup_target(
                    tables, target, charge, multiplicity
                )
            if target_entry is not None:
                return self._reference_species_result(
                    target, target_entry, tables, via=target_via
                )
            if resolver is not None:
                raise RuntimeError(
                    "No CBH rung has all its reference species in ATcT "
                    f"(ver. {tables.version}) or the reference library; see the "
                    "reactions above. Supply the missing ones as your own "
                    "references ('seamm-thermochemistry-references', with 'Use my "
                    "own references'), or compute them."
                )
            raise RuntimeError(
                "No CBH rung has all its reference species in ATcT "
                f"(ver. {tables.version}); see the reactions above."
            )

        # The species to compute: the target and the unique fragments of the
        # usable rungs, each once
        fragments = cbh.unique_fragments([reactions[r] for r in usable])

        mode = _P["cbh energies"]
        # The species are optimized over the warm engine, so the CBH approach
        # drives the Model Chemistry as an engine unless told otherwise
        route = _P.get("hessian route", "automatic") or "automatic"
        if route == "automatic":
            route = "warm engine over MDI"
        self._hessian_route = route
        if mode != "existing results (run nothing)":
            text = f"The Model Chemistry is driven as: {route}.\n"
            printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))
        energies = {}
        sources = {}
        configurations = {target: target_configuration}
        store = self._fragment_store(_P)
        energy_mc = self._energy_level(_P)
        if energy_mc is not None:
            if mode == "existing results (run nothing)":
                raise RuntimeError(
                    "Energies at another level than the geometries need the species "
                    "computed: choose 'finite-difference harmonic' for the species' "
                    "energies, or leave 'Energies at' as the geometry level."
                )
            if _P["thermal corrections"] != "harmonic":
                raise RuntimeError(
                    "Energies at another level than the geometries need harmonic "
                    "thermal corrections at the geometry level."
                )
            printer.important(
                __(
                    f"The energies are at {energy_mc['level']}, on the geometries "
                    f"and with the frequencies of {self.model}.\n",
                    indent=4 * " ",
                    wrap=False,
                    dedent=False,
                )
            )
        if mode == "existing results (run nothing)":
            energies[target] = self._energies_from_properties(
                _P, target_configuration, target
            )
            sources[target] = "properties"
            for can in fragments:
                configuration = self._find_fragment(system_db, can)
                if configuration is None:
                    stored = self._energies_from_store(_P, store, can)
                    if stored is None:
                        raise RuntimeError(
                            f"No system named 'CBH:{can}' with results for the "
                            f"fragment {can}, and nothing in the fragment store; "
                            "run the fragments first or use another mode."
                        )
                    energies[can], sources[can] = stored
                else:
                    energies[can] = self._energies_from_properties(
                        _P, configuration, can
                    )
                    sources[can] = "properties"
        else:
            energies[target] = self._compute_species(
                _P,
                target,
                smiles,
                target_system,
                target_configuration,
                embed=str(_P.get("target conformer", "as given")).startswith("lowest"),
            )
            sources[target] = "computed"
            for can, generated in fragments.items():
                # Already computed in this job (e.g. by an earlier iteration of a
                # Loop), with the fingerprints expected for it?
                in_job = self._from_job(_P, system_db, can, energy_mc is not None)
                if in_job is not None:
                    energies[can], configurations[can] = in_job
                    sources[can] = "job"
                    printer.important(
                        __(
                            f"\nSpecies {can}: already computed in this job.",
                            indent=4 * " ",
                            wrap=False,
                            dedent=False,
                        )
                    )
                    continue
                if energy_mc is not None:
                    found = self._geometry_from_store(store, can)
                    if found is not None:
                        energies[can], geometry, sources[can] = found
                        configurations[can] = self._configuration_from_geometry(
                            system_db, can, geometry
                        )
                        printer.important(
                            __(
                                f"\nSpecies {can}: geometry and frequencies from the "
                                f"fragment store ({sources[can]}).",
                                indent=4 * " ",
                                wrap=False,
                                dedent=False,
                            )
                        )
                        continue
                    energies[can] = self._compute_species(_P, can, generated)
                    sources[can] = "computed"
                    configurations[can] = self._find_fragment(system_db, can)
                    continue
                stored = self._energies_from_store(_P, store, can)
                if stored is not None:
                    energies[can], sources[can] = stored
                    how = "unverified" if "unverified" in sources[can] else "verified"
                    printer.important(
                        __(
                            f"\nSpecies {can}: from the fragment store ({how}).",
                            indent=4 * " ",
                            wrap=False,
                            dedent=False,
                        )
                    )
                    continue
                energies[can] = self._compute_species(_P, can, generated)
                sources[can] = "computed"
            # The energy level's single points on the geometries
            if energy_mc is not None:
                self._energy_level_energies(
                    energy_mc,
                    store,
                    target,
                    [target] + list(fragments),
                    configurations,
                    energies,
                    sources,
                )
            self._report_store_notes(store)
            # Add what was computed here to the fragment store (write-back)
            computed = {}
            for can in [target] + list(fragments):
                if "computed" not in sources.get(can, ""):
                    continue
                configuration = configurations.get(can)
                if configuration is None:
                    configuration = self._find_fragment(system_db, can)
                computed[can] = configuration
            self._write_back(_P, computed)
            # Back to the target
            system_db.system = target_system
            target_system.configuration = target_configuration
            self._working_configuration = target_configuration

        # The ZPE scale factor
        zpe_scale, zpe_source, zpe_flag = self._zpe_scale(_P)
        if zpe_scale != 1.0 or zpe_flag:
            text = f"The ZPEs are scaled by {zpe_scale} ({zpe_source})"
            text += f"; note: {zpe_flag}." if zpe_flag else "."
            printer.important(__(text, indent=4 * " ", wrap=True, dedent=False))
        if zpe_scale != 1.0:
            for e in energies.values():
                if "ZPE" not in e:
                    continue
                shift = (zpe_scale - 1.0) * e["ZPE"]
                e["ZPE scaled"] = zpe_scale * e["ZPE"]
                if "H0" in e:
                    e["H0"] += shift
                if "H298" in e:
                    e["H298"] += shift

        # Combine
        data = {}
        rows = []
        table = self.get_table("CBH")
        for column, coltype in (
            ("rung", "integer"),
            ("reaction", "string"),
            ("DfH(0 K)", "float"),
            ("DfH(298.15 K)", "float"),
            ("uncertainty (95%)", "float"),
            ("uncertainty method", "string"),
            ("DrH(0 K)", "float"),
            ("DrH(298.15 K)", "float"),
            ("ATcT DfH(0 K)", "float"),
            ("ATcT DfH(298.15 K)", "float"),
            ("ATcT uncertainty (95%)", "float"),
            ("ATcT version", "string"),
        ):
            if column not in table.columns:
                table.add_column(column, coltype)

        # A library value given at one temperature only is converted to the
        # other with this run's own H298 - H0 of the species (science: never a
        # cached or other-level value), and says so
        converted = {}

        def ref_value(can, column):
            entry = entries[can]
            value = entry.get(column)
            if value is not None:
                return value
            other = entry.get("dfH0_298K" if column == "dfH0_0K" else "dfH0_0K")
            e = energies.get(can, {})
            if other is None or "H298" not in e or "H0" not in e:
                return None
            if not entry.get("formula"):
                return None
            from seamm_thermochemistry.resolver import dfh_0k_from_298k

            thermal = e["H298"] - e["H0"]
            elements = dfh_0k_from_298k(0.0, 0.0, entry["formula"])
            if column == "dfH0_0K":
                value = other - thermal + elements
                converted[can] = ("0 K", value, thermal)
            else:
                value = other + thermal - elements
                converted[can] = ("298.15 K", value, thermal)
            return value

        results = {}
        fetch_covariances = _P["covariances"] == "fetch if needed"
        fetch_failed = False
        for rung in usable:
            reaction = reactions[rung]
            result = {"rung": rung, "reaction": str(reaction)}
            for key, label in (("H0", "0 K"), ("H298", "298.15 K")):
                if all(key in energies[s] for s in reaction.coefficients):
                    drH = sum(
                        c * energies[s][key] for s, c in reaction.coefficients.items()
                    )
                    column = "dfH0_0K" if key == "H0" else "dfH0_298K"
                    values = {s: ref_value(s, column) for s in reaction.references}
                    if any(v is None for v in values.values()):
                        continue
                    ref = sum(c * values[s] for s, c in reaction.references.items())
                    result[f"DrH {label}"] = drH
                    result[f"DfH {label}"] = ref - drH
            # The ATcT part: the reaction's coefficients on ATcT species,
            # through the ladder's coefficients where a reference is laddered
            atct_c = {}
            method_c2 = 0.0
            for s, c in reaction.references.items():
                for _id, cc in entries[s]["atct_coefficients"].items():
                    atct_c[_id] = atct_c.get(_id, 0.0) + c * cc
                if entries[s].get("source") == "laddered":
                    method_c2 += c * c
            atct_c = {k: v for k, v in atct_c.items() if abs(v) > 1e-12}
            ids = list(atct_c)
            coefficients = list(atct_c.values())
            uncertainties = [
                tables.entries[_id]["uncertainty_298K"] or 0.0 for _id in ids
            ]
            # Fetch a missing block at most until a fetch fails (each failure
            # waits for the timeout, e.g. on a compute node without network)
            fetch = fetch_covariances and not fetch_failed
            covariance = tables.covariance(ids, fetch=fetch) if ids else None
            if ids and covariance is None and fetch:
                fetch_failed = True
            u_atct, used = reaction_uncertainty(coefficients, uncertainties, covariance)
            u_method = 0.0
            if method_c2 > 0 and ladder is not None:
                rms = ladder.method_uncertainty.get("rms_0K")
                if rms:
                    u_method = rms * method_c2**0.5
            # Library references: no covariance with ATcT, in quadrature
            u_library = 0.0
            for s, c in reaction.references.items():
                if entries[s].get("source") not in ("ATcT", "laddered"):
                    u_library += (c * (entries[s]["uncertainty_298K"] or 0.0)) ** 2
            u_library = u_library**0.5
            u = (u_atct**2 + u_method**2 + u_library**2) ** 0.5
            result["uncertainty"] = u
            result["uncertainty ATcT"] = u_atct
            if used:
                result["uncertainty method"] = "ATcT covariance"
            elif fetch_covariances:
                result["uncertainty method"] = "quadrature (covariance fetch failed)"
            else:
                result["uncertainty method"] = "quadrature (no covariance cached)"
            if u_method:
                result["uncertainty ladder method"] = u_method
                result["uncertainty method"] += " + ladder method term"
            if u_library:
                result["uncertainty reference library"] = u_library
                result["uncertainty method"] += " + reference library (quadrature)"
            results[rung] = result
            rows.append(result)
            table.append_row(
                **{
                    "rung": rung,
                    "reaction": str(reaction),
                    "DfH(0 K)": result.get("DfH 0 K"),
                    "DfH(298.15 K)": result.get("DfH 298.15 K"),
                    "uncertainty (95%)": u,
                    "uncertainty method": result["uncertainty method"],
                    "DrH(0 K)": result.get("DrH 0 K"),
                    "DrH(298.15 K)": result.get("DrH 298.15 K"),
                    "ATcT DfH(0 K)": target_entry["dfH0_0K"] if target_entry else None,
                    "ATcT DfH(298.15 K)": (
                        target_entry["dfH0_298K"] if target_entry else None
                    ),
                    "ATcT uncertainty (95%)": (
                        target_entry["uncertainty_298K"] if target_entry else None
                    ),
                    "ATcT version": tables.version,
                }
            )

        # Print the species and the results
        species_rows = []
        for can in [target] + list(fragments):
            e = energies[can]
            entry = entries.get(can) if can != target else target_entry
            species_rows.append(
                [
                    can if can != target else f"{can} (target)",
                    sources.get(can, ""),
                    f"{e['E']:.3f}" if "E" in e else "",
                    (f"{e['E geometry level']:.3f}" if "E geometry level" in e else ""),
                    f"{e['ZPE']:.3f}" if "ZPE" in e else "",
                    f"{e['H298'] - e['E']:.3f}" if "E" in e and "H298" in e else "",
                    (
                        (
                            f"ladder (CBH-{entry['rung']})"
                            if entry.get("source") == "laddered"
                            else (
                                entry["id"]
                                if entry.get("source", "ATcT") == "ATcT"
                                else entry["source"]
                            )
                        )
                        if entry
                        else ""
                    ),
                    _reference_cell(entry, converted.get(can), "0 K", "dfH0_0K"),
                    _reference_cell(entry, converted.get(can), "298.15 K", "dfH0_298K"),
                    (
                        "exact"
                        if entry and entry["exact"]
                        else f"{entry['uncertainty_298K']:.3f}" if entry else ""
                    ),
                ]
            )
        used_library = {
            can: entries[can]
            for rung in usable
            for can in reactions[rung].references
            if entries[can].get("source") not in ("ATcT", "laddered")
        }
        ref_name = "Ref." if used_library else "ATcT"
        text = "\nThe species (energies in kJ/mol):\n"
        text += tabulate(
            species_rows,
            headers=[
                "Species",
                "Source",
                "E",
                "E (geometry level)",
                "ZPE",
                "H(298.15 K) - E",
                "Reference",
                f"{ref_name} ΔfH°(0 K)",
                f"{ref_name} ΔfH°(298.15 K)",
                "± (95%)",
            ],
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=(
                "left",
                "left",
                "decimal",
                "decimal",
                "decimal",
                "decimal",
                "left",
                "decimal",
                "decimal",
                "decimal",
            ),
        )
        text += "\n"

        result_rows = []
        for rung in usable:
            r = results[rung]
            row = [f"CBH-{rung}"]
            for key in ("DrH 0 K", "DfH 0 K", "DrH 298.15 K", "DfH 298.15 K"):
                row.append(f"{r[key]:.2f}" if key in r else "")
            row.append(f"{r['uncertainty']:.2f}")
            if target_entry is not None:
                for key, column in (
                    ("DfH 0 K", "dfH0_0K"),
                    ("DfH 298.15 K", "dfH0_298K"),
                ):
                    if key in r and target_entry[column] is not None:
                        row.append(f"{r[key] - target_entry[column]:+.2f}")
                    else:
                        row.append("")
            result_rows.append(row)
        headers = [
            "Rung",
            "ΔrH(0 K)",
            "ΔfH°(0 K)",
            "ΔrH(298.15 K)",
            "ΔfH°(298.15 K)",
            "± ref (95%)",
        ]
        if target_entry is not None:
            headers += ["Δ vs ATcT (0 K)", "Δ vs ATcT (298.15 K)"]
        text += "\nEnthalpy of formation of the target (kJ/mol):\n"
        text += tabulate(
            result_rows,
            headers=headers,
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left",) + ("decimal",) * (len(headers) - 1),
        )
        text += "\n"
        if target_entry is not None:
            unc = (
                "exact"
                if target_entry["exact"]
                else f"± {target_entry['uncertainty_298K']:.2f}"
            )
            text += (
                f"\nATcT ({target_entry['id']}, ver. {tables.version}) gives "
                f"ΔfH°(0 K) = {target_entry['dfH0_0K']:.2f} and ΔfH°(298.15 K) = "
                f"{target_entry['dfH0_298K']:.2f} {unc} kJ/mol for the target.\n"
            )
        if used_library:
            text += "\nReference values from the reference library:\n"
            for can, entry in used_library.items():
                text += f"    {can}: {entry['source']}, {entry['citation']}"
                if entry.get("override_reason"):
                    text += f" (replaces an earlier source: {entry['override_reason']})"
                text += ".\n"
            if converted:
                text += (
                    "* Converted between 0 K and 298.15 K with this calculation's "
                    "own H(298.15 K) - H(0) of the species (0 K via computed "
                    "H298-H0).\n"
                )
            text += (
                "Library values have no covariance with ATcT; their uncertainties "
                "(95 %, with each source's offset against ATcT) are added in "
                "quadrature. The source offsets come from well-studied species and "
                "are a lower bound for these.\n"
            )
        text += (
            "\nThe reference uncertainty is the 95 % confidence limit of the ATcT "
            "part only; ATcT uncertainties and covariances are defined at 298.15 K "
            "and are used for the 0 K result as an approximation. The computed "
            "reaction enthalpy adds its own, unquantified, error. A plain ATcT "
            "gas-phase entry is the lowest conformer at 0 K but a conformer average "
            "at 298.15 K, which limits the 298.15 K comparison for flexible molecules."
        )
        if any(src.startswith("store (unverified") for src in sources.values()):
            text += (
                "\nSome species came from the fragment store without a settings "
                "fingerprint (matched on the model chemistry only), so they are "
                "unverified: check that their settings match this calculation."
            )
        methods = [r["uncertainty method"] for r in results.values()]
        if any("fetch failed" in m for m in methods):
            text += (
                "\nFetching the ATcT covariance blocks from atct.anl.gov failed "
                "(see the log), so some reference uncertainties ignore the "
                "correlations (quadrature); run 'seamm-thermochemistry-atct "
                "covariance' on a machine with network access to cache them."
            )
        elif any(m.startswith("quadrature") for m in methods):
            text += (
                "\nNo ATcT covariance block was cached for some reactions, so their "
                "reference uncertainties ignore the correlations (quadrature); run "
                "'seamm-thermochemistry-atct covariance' to cache them, or set "
                "'ATcT covariances' to fetch them."
            )
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        # The spread between the rungs (CBH-1 up), an indicator of the uncertainty
        spread, spread_rungs = rung_spread(results, usable)
        if spread is not None:
            printer.important(
                __(
                    "\nThe spread of ΔfH°(298.15 K) over CBH-"
                    + ", CBH-".join(str(r) for r in spread_rungs)
                    + f" is {spread:.2f} kJ/mol. A large spread points to trouble "
                    "with the method or the reference values.",
                    indent=4 * " ",
                    wrap=True,
                    dedent=False,
                )
            )

        # The results: the highest usable rung
        best = results[usable[-1]]
        if preset == "quick" and best["rung"] < 2:
            printer.important(
                __(
                    f"\nWarning: the highest rung available is CBH-{best['rung']}, "
                    "where r2SCAN-3c is weak (several kJ/mol at CBH-1; aromatic "
                    "rings and small rings are kept whole, so small aromatics stop "
                    "below CBH-2). Use the 'default' preset for this molecule.",
                    indent=4 * " ",
                    wrap=True,
                    dedent=False,
                )
            )
        data["ZPE scale factor"] = zpe_scale
        data["ZPE scale factor source"] = zpe_source + (
            f" ({zpe_flag})" if zpe_flag else ""
        )
        if spread is not None:
            data["CBH rung spread"] = spread
        data["CBH energy level"] = (
            energy_mc["level"] if energy_mc is not None else self.model
        )
        data["CBH geometry level"] = self.model
        composite = getattr(self, "_composite", None)
        if composite is not None:
            data["CBH composite recipe"] = composite["recipe"]
            data["CBH composite recipe version"] = composite["version"]
            data["CBH composite recipe fingerprint"] = composite["fingerprint"]
        if getattr(self, "_high_choice", None) is not None:
            data["CBH high level"] = self._high_choice["reason"]
            data["CBH high canonical limit"] = int(self._high_max_bf)
        data["CBH rung"] = best["rung"]
        data["CBH reaction"] = best["reaction"]
        if "DfH 0 K" in best:
            data["DfH 0K"] = best["DfH 0 K"]
        if "DfH 298.15 K" in best:
            data["DfH 298K"] = best["DfH 298.15 K"]
        data["DfH uncertainty"] = best["uncertainty"]
        data["DfH uncertainty method"] = best["uncertainty method"]
        data["ATcT version"] = tables.version
        if resolver is not None:
            best_refs = {
                can: entries[can] for can in reactions[best["rung"]].references
            }
            from seamm_thermochemistry.references import elements

            several = {}
            for can, e in best_refs.items():
                if e.get("source") in ("ATcT", "laddered"):
                    continue
                for symbol in elements(can):
                    if not resolver.atct_has_element(symbol):
                        several.setdefault(symbol, set()).add(e["source"])
            several = {k: sorted(v) for k, v in several.items() if len(v) > 1}
            data["CBH references"] = resolver.label(best_refs, several)
            data["references"] = {
                **resolver.versions(),
                "library values": {
                    can: {
                        "id": e["id"],
                        "source": e["source"],
                        "citation": e["citation"],
                        "converted": (
                            f"{converted[can][0]} via computed H298-H0"
                            if can in converted
                            else None
                        ),
                    }
                    for can, e in best_refs.items()
                    if e.get("source") not in ("ATcT", "laddered")
                },
            }
        else:
            data["CBH references"] = "ATcT-only" if ladder is None else "laddered"
        if target_entry is not None:
            data["ATcT DfH 0K"] = target_entry["dfH0_0K"]
            data["ATcT DfH 298K"] = target_entry["dfH0_298K"]
            data["ATcT DfH uncertainty"] = target_entry["uncertainty_298K"]
        data["CBH results"] = {
            str(rung): {k: v for k, v in r.items()} for rung, r in results.items()
        }

        # Citations
        self._cite_cbh(tables)
        return data

    def _ensure_atct(self, base, version):
        """Download ATcT when no version is installed (e.g. an installation
        made before the installer set it up). Returns True if it downloaded.

        atct.anl.gov serves only its current version, so a specific version
        is never downloaded. Parallel jobs may race to download:
        the download goes to a per-process file and an existing version wins."""
        from seamm_thermochemistry import atct

        if atct.ATcT.available(base):
            return False
        directory = atct.atct_directory(base)
        if version is not None:
            raise RuntimeError(
                f"ATcT version {version} is not installed in {directory}. Install "
                "it with 'seamm-thermochemistry-atct download --from FILE' from a "
                "copy of that version, or choose 'latest'."
            )
        try:
            atct.download_all(base, exist_ok=True)
        except (atct.ATcTError, OSError) as e:
            raise RuntimeError(
                f"No ATcT data are installed in {directory}, and downloading them "
                f"failed: {e}. Install them with 'seamm-thermochemistry-installer "
                "install' (or 'seamm-thermochemistry-atct download --from FILE' "
                "from a copy saved in a browser) on this machine, or set the ATcT "
                "directory to an installed copy."
            ) from None
        return True

    def _reference_species_result(self, target, entry, tables, via=None):
        """The result for a target that is itself a reference species: its ATcT
        value, with nothing computed. `via`: the ATcT SMILES it matched, when
        ATcT writes another resonance form."""
        dfh0, dfh298 = entry["dfH0_0K"], entry["dfH0_298K"]
        u = 0.0 if entry["exact"] else (entry["uncertainty_298K"] or 0.0)
        label = f"target is a reference species (ATcT {entry['id']})"
        text = (
            f"\nNo rung gives a reaction for {target}, but the target is itself a "
            f"reference species, ATcT {entry['id']} (ver. {tables.version}"
        )
        if via is not None:
            text += (
                f", written {via}: the same species as {target}, another resonance "
                "form"
            )
        text += "): "
        if dfh0 is not None:
            text += f"ΔfH°(0 K) = {dfh0:.2f} and "
        text += f"ΔfH°(298.15 K) = {dfh298:.2f} "
        text += "(exact)" if entry["exact"] else f"± {u:.2f}"
        text += " kJ/mol. That value is the result; no species were computed."
        printer.important(__(text, indent=4 * " ", wrap=True, dedent=False))
        data = {
            "CBH reaction": label,
            "DfH 298K": dfh298,
            "DfH uncertainty": u,
            "DfH uncertainty method": "the target's own ATcT value",
            "ATcT version": tables.version,
            "ATcT DfH 298K": dfh298,
            "ATcT DfH uncertainty": entry["uncertainty_298K"],
        }
        if dfh0 is not None:
            data["DfH 0K"] = dfh0
            data["ATcT DfH 0K"] = dfh0
        self._cite_cbh(tables)
        return data

    def _run_species_only(self, _P, smiles, target, system, configuration):
        """Compute only the target, exactly as a CBH species: optimized and its
        Hessian at the geometry level, then the energy level's single point,
        with the properties and fingerprints the fragment store's harvest reads.
        No reactions and no ATcT. Returns the data dictionary."""
        preset = _P.get("preset", "custom")
        if _P["cbh energies"] == "existing results (run nothing)":
            raise RuntimeError(
                "'The species only' computes the species: choose "
                "'finite-difference harmonic' for the species' energies."
            )
        text = ""
        if preset in PRESETS:
            levels = PRESETS[preset]
            text += (
                f"The '{preset}' preset ({levels['recipe']}): the geometry and "
                "frequencies with "
                f"{levels['geometry']}, the energy with "
                f"{levels['energy'] or levels['geometry']}.\n"
            )
        text += f"Only the species {smiles} is computed, for the fragment store"
        charge = configuration.charge
        multiplicity = configuration.spin_multiplicity or 1
        if charge != 0 or multiplicity != 1:
            text += f" (charge {charge}, multiplicity {multiplicity})"
        text += ": no CBH reactions and no ATcT.\n"
        route = _P.get("hessian route", "automatic") or "automatic"
        if route == "automatic":
            route = "warm engine over MDI"
        self._hessian_route = route
        text += f"The Model Chemistry is driven as: {route}.\n"
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        store = self._fragment_store(_P)
        energy_mc = self._energy_level(_P)
        if energy_mc is not None and _P["thermal corrections"] != "harmonic":
            raise RuntimeError(
                "Energies at another level than the geometries need harmonic "
                "thermal corrections at the geometry level."
            )
        energies = {
            target: self._compute_species(
                _P,
                target,
                smiles,
                system,
                configuration,
                embed=str(_P.get("target conformer", "as given")).startswith("lowest"),
            )
        }
        sources = {target: "computed"}
        if energy_mc is not None:
            printer.important(
                __(
                    f"\nThe energy is at {energy_mc['level']}, on the geometry and "
                    f"with the frequencies of {self.model}.",
                    indent=4 * " ",
                    wrap=False,
                    dedent=False,
                )
            )
            self._energy_level_energies(
                energy_mc,
                store,
                target,
                [target],
                {target: configuration},
                energies,
                sources,
            )
        self._report_store_notes(store)
        self._write_back(_P, {target: configuration})
        system_db = self.get_variable("_system_db")
        system_db.system = system
        system.configuration = configuration
        self._working_configuration = configuration

        e = energies[target]
        row = [target, sources[target], f"{e['E']:.3f}" if "E" in e else ""]
        row.append(f"{e['E geometry level']:.3f}" if "E geometry level" in e else "")
        row.append(f"{e['ZPE']:.3f}" if "ZPE" in e else "")
        row.append(f"{e['H298'] - e['E']:.3f}" if "E" in e and "H298" in e else "")
        text = "\nThe species (energies in kJ/mol, ZPE unscaled):\n"
        text += tabulate(
            [row],
            headers=[
                "Species",
                "Source",
                "E",
                "E (geometry level)",
                "ZPE",
                "H(298.15 K) - E",
            ],
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "left", "decimal", "decimal", "decimal", "decimal"),
        )
        if e.get("imaginary"):
            text += (
                f"\nWarning: the structure has {e['imaginary']} imaginary mode(s): "
                "it is not a minimum, and the store's lookups leave it out."
            )
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        data = {
            "CBH energy level": (
                energy_mc["level"] if energy_mc is not None else self.model
            ),
            "CBH geometry level": self.model,
        }
        composite = getattr(self, "_composite", None)
        if composite is not None:
            data["CBH composite recipe"] = composite["recipe"]
            data["CBH composite recipe version"] = composite["version"]
            data["CBH composite recipe fingerprint"] = composite["fingerprint"]
        return data

    # -----------------------------------------------------------------------

    def _embed_lowest(self, configuration, smiles, n=20):
        """Rebuild the configuration's atoms and coordinates as the lowest
        MMFF94 conformer of `n` RDKit embeddings of `smiles`, keeping its
        charge and spin multiplicity, and say what was found."""
        from seamm_thermochemistry import cbh

        mol, symbols, xyz, energies, best = cbh.lowest_conformer(
            smiles, n=n, return_mol=True
        )
        charge = configuration.charge
        multiplicity = configuration.spin_multiplicity
        configuration.from_RDKMol(mol, properties="none")
        configuration.charge = charge
        configuration.spin_multiplicity = multiplicity
        spread = max(energies) - min(energies)
        printer.important(
            __(
                f"    starting from the lowest of {len(energies)} MMFF94 conformers "
                f"(spread {spread:.1f} kcal/mol)" + self._chain_dihedrals(symbols, xyz),
                indent=4 * " ",
                wrap=False,
                dedent=False,
            )
        )

    @staticmethod
    def _chain_dihedrals(symbols, xyz):
        """', heavy-atom chain dihedrals ...' for an unbranched chain of 4+
        heavy atoms, else ''."""
        import numpy as np

        heavy = [i for i, s in enumerate(symbols) if s != "H"]
        if len(heavy) < 4:
            return ""
        p = np.asarray(xyz, dtype=float)
        text = []
        for i in range(len(heavy) - 3):
            a, b, c, d = (p[heavy[i + k]] for k in range(4))
            b0, b1, b2 = a - b, c - b, d - c
            b1n = b1 / np.linalg.norm(b1)
            v = b0 - np.dot(b0, b1n) * b1n
            w = b2 - np.dot(b2, b1n) * b1n
            angle = np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))
            text.append(f"{abs(angle):.0f}")
        return ", heavy-atom dihedrals (in SMILES order) " + " ".join(text) + "°"

    #: The "fragment store" choices that read the store ("use if available" is
    #: the name of "use only" in flowcharts from before 2026.11.12)
    _STORE_READ = ("use and add what is computed", "use only", "use if available")

    def _job_root(self):
        """The job's directory, or None outside a flowchart."""
        try:
            return Path(self.flowchart.root_directory)
        except Exception:
            return None

    def _job_label(self):
        """The job as the fragment store records it: '<job id>:<directory>'."""
        root = self._job_root()
        if root is None:
            return None
        job_id = None
        try:
            text = (root / "job_data.json").read_text()
            if text.startswith("!"):  # SEAMM's header, "!MolSSI job_data 1.0"
                text = text.split("\n", 1)[1] if "\n" in text else ""
            job_id = json.loads(text).get("job id")
        except Exception:
            pass
        return f"{job_id}:{root}" if job_id is not None else str(root)

    def _write_back(self, _P, species):
        """Add the species computed here ({SMILES: configuration}) to the
        fragment store (write-back), with the guards of a running job's
        results: a structure whose bonds are not its species', or with an
        imaginary mode, is refused; a lower, distinct conformer is added as one.

        Only with "use and add what is computed", and only when the store is on
        a local disk of this machine; otherwise the job is marked to be
        harvested when it comes back to its JobServer (or by hand with
        ``seamm-thermochemistry-fragments harvest``). Never fails the step."""
        if _P.get("fragment store", "ignore") != "use and add what is computed":
            return
        species = {k: v for k, v in species.items() if v is not None}
        if not species:
            return
        try:
            from seamm_thermochemistry import (
                FragmentStore,
                harvest_configuration,
                hooks,
                store_is_writable,
            )
            from seamm_thermochemistry.fragment_store import new_counts
        except ImportError:  # pragma: no cover - an older seamm-thermochemistry
            return
        path = _P.get("fragment store file", "default")
        path = None if path in ("default", "") else os.path.expanduser(path)
        root = self._job_root()

        def not_written(why):
            if root is not None:
                hooks.write_marker(root, list(species), why)
            text = (
                f"\nThe computed species were not added to the fragment store: {why}. "
                "The job will be harvested into the store when it returns to its "
                "JobServer, or harvest it with 'seamm-thermochemistry-fragments "
                "harvest'."
            )
            printer.important(__(text, indent=4 * " ", wrap=True, dedent=False))

        writable, why = store_is_writable(path)
        if not writable:
            not_written(f"the store is {why}")
            return
        counts = new_counts()
        job = self._job_label()
        try:
            with FragmentStore(path, create=True) as store:
                for can, configuration in species.items():
                    harvest_configuration(
                        store,
                        configuration,
                        can,
                        charge=cbh_charge(can),
                        multiplicity=cbh_multiplicity(can),
                        job=job,
                        energy="total energy#Thermochemistry#*",
                        zpe="ZPE#Thermochemistry#*",
                        enthalpy="H(298.15 K)#Thermochemistry#*",
                        fingerprint="settings fingerprint#Thermochemistry#*",
                        frequency_fingerprint=(
                            "frequency fingerprint#Thermochemistry#*"
                        ),
                        write_back=True,
                        counts=counts,
                    )
                where = store.path
        except Exception as e:
            logger.warning(f"Writing the fragment store failed: {e}")
            not_written(f"writing it failed ({e})")
            return

        def tally(table):
            return ", ".join(f"{n} {k}" for k, n in sorted(counts[table].items()))

        text = f"\nAdded to the fragment store ({where}): "
        text += f"geometries: {tally('geometry') or 'none'}; "
        text += f"energies: {tally('energy') or 'none'}."
        for reason in counts["skipped"]:
            if "no property matching" in reason:
                continue  # a species with nothing computed at a level
            text += f"\n    {reason}"
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

    def _fragment_store(self, _P):
        """The fragment store to read, or None if ignored or unreadable."""
        if _P.get("fragment store", "ignore") not in self._STORE_READ:
            return None
        path = _P.get("fragment store file", "default")
        path = None if path in ("default", "") else os.path.expanduser(path)
        self._store_layer = {}
        try:
            # The local store, then the official one (seamm-thermochemistry's
            # two layers); an older seamm-thermochemistry has only the local one
            try:
                from seamm_thermochemistry import LayeredStore
            except ImportError:  # pragma: no cover
                from seamm_thermochemistry import FragmentStore

                return FragmentStore(path)
            return LayeredStore(path)
        except Exception as e:  # pragma: no cover
            logger.warning(f"The fragment store could not be opened: {e}")
            return None

    @staticmethod
    def _store_source(store, row):
        """How a value from the store is labelled: "store", with the layer when
        it is the official one, and "unverified" when it has no fingerprint."""
        layer = row.get("layer", "local")
        text = "store"
        details = []
        if layer != "local" and hasattr(store, "layer_name"):
            details.append(store.layer_name(layer))
        if not row.get("verified", True):
            details.append("unverified")
        return text + (f" ({', '.join(details)})" if details else "")

    def _report_store_notes(self, store):
        """Print what the two-layer store noticed (an official value superseding
        an old local one, local and official values that differ)."""
        notes = getattr(store, "notes", None)
        if not notes:
            return
        text = "\nFragment store: " + "\n    ".join(notes)
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))
        notes.clear()

    def _energies_from_store(self, _P, store, can):
        """The energies of a species from the store at this model chemistry,
        as ({"E", "ZPE", "H0", "H298"}, source), or None. A record must have
        what the thermal-correction mode needs: E and H(298.15 K) with
        harmonic corrections (ZPE too for the 0 K result), E alone without."""
        if store is None:
            return None
        from thermochemistry_step.cbh_approach import cbh_charge, cbh_multiplicity

        model = self.model or "unknown"
        # The fingerprints expected for this species, not the target's
        energy_fingerprints, _ = self._geometry_level_fingerprints(can)
        record = store.get(
            can,
            charge=cbh_charge(can),
            multiplicity=cbh_multiplicity(can),
            energy_fingerprint=energy_fingerprints,
            model_chemistry=model,
        )
        if record is None or record.get("E") is None:
            return None
        result = {"E": record["E"]}
        if _P["thermal corrections"] != "harmonic":
            result["H298"] = record["E"]
        else:
            if record.get("H") is None:
                return None
            if record.get("T") not in (None, T_REF):
                return None
            result["H298"] = record["H"]
            if record.get("ZPE") is not None:
                result["ZPE"] = record["ZPE"]
                result["H0"] = record["E"] + record["ZPE"]
        return result, self._store_source(store, record)

    @staticmethod
    def _preset_rungs(_P):
        """The rungs a preset uses: up to CBH-2, or CBH-3 too."""
        rungs = [0, 1, 2]
        if _yes(_P.get("include CBH-3", "no")):
            rungs.append(3)
        return rungs

    def apply_preset(self, _P):
        """The parameters with the preset's settings in place: the rungs, harmonic
        thermal corrections, computed species and the energy level. 'custom'
        changes nothing. The geometry level is the preset's model chemistry
        (see :meth:`preset_model_chemistry`)."""
        preset = _P.get("preset", "custom")
        if preset not in PRESETS:
            return _P
        P = dict(_P)
        P["cbh rungs"] = ", ".join(str(r) for r in self._preset_rungs(_P))
        P["thermal corrections"] = "harmonic"
        P["cbh energies"] = "finite-difference harmonic"
        energy = self.preset_levels(preset)["energy"]
        P["energy model chemistry"] = energy or "the geometry level"
        return P

    def preset_levels(self, preset):
        """The geometry and energy levels of a preset: its recipe's, when the
        Composite step is installed, else the spellings in :data:`PRESETS`."""
        levels = dict(PRESETS[preset])
        if levels.get("recipe") is None:
            return levels
        try:
            import composite_step  # noqa: F401
        except ImportError:
            return levels
        recipe = recipe_levels(levels["recipe"], getattr(self, "data_path", None))
        levels.update(
            geometry=recipe["geometry"], energy=recipe["energy"], composite=recipe
        )
        return levels

    def preset_model_chemistry(self, _P):
        """The ``_model_chemistry`` wrapper of the preset's geometry level, or None
        for 'custom'. A preset's recipe (or a composite recipe as the flowchart's
        model chemistry) is kept in ``self._composite``: its geometry level is the
        step's model, its energy level and ZPE factor follow."""
        preset = _P.get("preset", "custom")
        self._composite = None
        if preset not in PRESETS:
            return self.composite_model_chemistry()
        levels = self.preset_levels(preset)
        self._composite = levels.get("composite")
        return self._resolve_model_chemistry(levels["geometry"])

    def composite_model_chemistry(self):
        """For a composite recipe as the flowchart's model chemistry, the wrapper
        of its geometry level (kept in ``self._composite``); else None."""
        import seamm

        if seamm.flowchart_variables is None:  # no flowchart (the GUI, tests)
            return None
        if not self.variable_exists("_model_chemistry"):
            return None
        level = (self.get_variable("_model_chemistry") or {}).get("level")
        if not is_composite(level):
            return None
        self._composite = recipe_levels(level, getattr(self, "data_path", None))
        return self._resolve_model_chemistry(self._composite["geometry"])

    def _zpe_scale(self, _P):
        """(factor, source, caution) for the ZPEs: a number given, or the table's
        factor for the geometry level, or 1 with a caution if it has none."""
        value = str(_P.get("ZPE scale factor", "from the table") or "").strip()
        if value and value != "from the table":
            try:
                return float(value), "as given", ""
            except ValueError:
                raise ValueError(
                    f"The ZPE scale factor must be 'from the table' or a number, "
                    f"not '{value}'."
                )
        from seamm_thermochemistry import zpe_scale_factor

        factor = zpe_scale_factor(self.model)
        composite = getattr(self, "_composite", None)
        if composite is not None and composite["geometry"] == self.model:
            # The recipe pins its factor; the same as the table's, it is reported
            # as the table's, so the results are those of the levels themselves
            scale, source = composite["zpe"]
            if factor is None or scale != factor.factor:
                return scale, f"the recipe {composite['name']}: {source}", ""
        if factor is None:
            return (
                1.0,
                "no factor in the table",
                f"the table has no ZPE scale factor for {self.model}, so the ZPEs "
                "are not scaled",
            )
        return factor.factor, factor.source, factor.flag

    def _energy_level(self, _P):
        """The ``_model_chemistry`` wrapper of the energy level, or None when the
        energies are at the geometry level."""
        text = str(_P.get("energy model chemistry", "") or "").strip()
        if text == HIGH_ENERGY:
            text = self._high_energy_level()
        composite = getattr(self, "_composite", None)
        if text in ("", "the geometry level") and composite is not None:
            # A composite as the flowchart's model chemistry: its energy level
            text = composite["energy"] or ""
        self._energy_recipe = None
        if is_composite(text):
            # A composite as the energy level: its energy, on this step's
            # geometries
            recipe = recipe_levels(text, getattr(self, "data_path", None))
            if recipe["geometry"] != self.model:
                printer.important(
                    __(
                        f"Note: the composite {text} is defined with geometries at "
                        f"{recipe['geometry']}; its energy, {recipe['energy']}, is "
                        f"used here on the geometries of {self.model}.",
                        indent=4 * " ",
                    )
                )
            if recipe["energy"] == recipe["recipe"]:
                # Several calculations, combined: the composite runs as itself
                self._energy_recipe = recipe["object"]
                return self._resolve_model_chemistry(
                    recipe["recipe"], getattr(self, "data_path", None)
                )
            text = recipe["energy"] or recipe["geometry"]
        if text in ("", "the geometry level") or text == self.model:
            return None
        return self._resolve_model_chemistry(text)

    def _high_energy_level(self):
        """The 'high' preset's energy level for the target (high_level.py); the
        default preset's when there is no 'high' result, with the reason printed
        and kept for the results."""
        from .high_level import high_level

        smiles, multiplicity, atomic_numbers = self._high_target
        max_bf = int(getattr(self, "_high_max_bf", None) or CANONICAL_MAX_BF_DEFAULT)
        level, reason = high_level(smiles, multiplicity, atomic_numbers, max_bf=max_bf)
        if level is None:
            level = PRESETS["default"]["energy"]
            note = (
                f"No 'high' result for {smiles}: {reason}. The energies are at the "
                f"default preset's level, {level}."
            )
        else:
            note = f"The 'high' energy level for {smiles}: {level}, since {reason}."
        self._high_choice = {"level": level, "reason": note}
        printer.important(__(note, indent=4 * " "))
        return level

    @staticmethod
    def _resolve_model_chemistry(text, data_path=None):
        """The ``_model_chemistry`` wrapper for a model chemistry string.
        ``data_path`` (the node's) lets discovery see a user's or a job's
        composite recipes."""
        import seamm

        try:
            from model_chemistry_step.model_chemistry import (
                availability_problem,
                discover_model_chemistries,
                match_model_chemistry,
                resolve_level,
            )
        except ImportError:  # pragma: no cover
            raise RuntimeError(
                "The presets and energies at another model chemistry need the "
                "Model Chemistry step (model-chemistry-step) installed."
            )
        text = resolve_level(text, dict(seamm.flowchart_variables._data))
        available = discover_model_chemistries(periodic_only=False, data_path=data_path)
        mc = available.get(text) or match_model_chemistry(text, available)
        if mc is None:
            raise ValueError(availability_problem(text, available, periodic=False))
        return mc

    def _geometry_from_store(self, store, can):
        """A species' geometry, ZPE and H(298.15 K) at the geometry level from the
        fragment store's geometry table (the lowest conformer at that level), as
        ({"E", "ZPE", "H0", "H298"}, geometry row, source), or None."""
        if store is None:
            return None
        # The fingerprints expected for this species, not the target's
        _, frequency_fingerprints = self._geometry_level_fingerprints(can)
        rows = store.geometries(
            can,
            charge=cbh_charge(can),
            multiplicity=cbh_multiplicity(can),
            frequency_fingerprint=frequency_fingerprints,
            model_chemistry=self.model,
        )
        rows = [
            r
            for r in rows
            if r["E"] is not None
            and r["ZPE"] is not None
            and r["H"] is not None
            and r.get("T") in (None, T_REF)
        ]
        if not rows:
            return None
        # The lowest conformer of both layers, a tie going to the local one
        lowest = getattr(store, "lowest", None)
        g = lowest(rows) if lowest else min(rows, key=lambda r: r["E"] + r["ZPE"])
        result = {"E": g["E"], "ZPE": g["ZPE"], "H0": g["E"] + g["ZPE"], "H298": g["H"]}
        if not hasattr(self, "_store_layer"):
            self._store_layer = {}
        self._store_layer[can] = g.get("layer", "local")
        return result, g, self._store_source(store, g)

    def _configuration_from_geometry(self, system_db, can, geometry):
        """The system CBH:<can> with the stored geometry as a new configuration."""
        system = None
        for candidate in system_db.get_systems():
            if candidate.name == f"CBH:{can}":
                system = candidate
                break
        if system is None:
            system = system_db.create_system(name=f"CBH:{can}", make_current=True)
        configuration = system.create_configuration(
            name="from the fragment store", make_current=True
        )
        xyz = np.asarray(geometry["coordinates"], dtype=float)
        configuration.atoms.append(
            symbol=geometry["symbols"], x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2]
        )
        configuration.charge = cbh_charge(can)
        configuration.spin_multiplicity = cbh_multiplicity(can)
        return configuration

    def _energy_level_energies(
        self, mc, store, target, order, configurations, energies, sources
    ):
        """The energy level's single points on the species' geometries.

        The target's runs first, which gives the energy level's settings
        fingerprint; then the fragments' energies come from the fragment store's
        energy table where it has them on the same geometry and level, and the
        rest run as one batch. Each species' energies become E(energy level) and
        [H - E](geometry level); the energy, its settings fingerprint and the
        calculation record are stored as properties on the configuration for the
        store's harvest."""
        import seamm_exec
        from seamm_thermochemistry.geometry_store import hash_geometry

        level = mc["level"]
        if is_composite(level):
            return self._composite_energies(
                mc, store, order, configurations, energies, sources
            )
        directory = self._species_directory()
        directory = directory / "energy_level" if directory is not None else None
        keys = {can: f"species_{i}" for i, can in enumerate(order)}
        by_key = {v: k for k, v in keys.items()}

        def run(cans):
            if not cans:
                return {}
            out = {}
            with seamm_exec.Evaluator(
                self,
                mc,
                properties=("energy",),
                directory=directory,
                name="Thermochemistry_energies",
            ) as evaluator:
                for can in cans:
                    evaluator.submit(configurations[can], key=keys[can])
                for result in evaluator.results():
                    out[by_key[result.key]] = result
            for can, result in out.items():
                if not result.ok or result.energy is None:
                    raise RuntimeError(
                        f"The {level} energy of {can} failed: {result.reason}"
                    )
            return out

        def record_of(result):
            record = (result.data or {}).get("calculation")
            return record if isinstance(record, dict) else None

        high = {}
        computed_parts = {}
        fingerprint = version = None
        for can, result in run([target]).items():
            high[can] = (result.energy, record_of(result), "computed")
            computed_parts[can] = (result.data or {}).get("energy parts")
            record = record_of(result)
            fingerprint = record.get("fingerprint") if record else None
            version = record.get("version") if record else None

        def expected(can):
            """The energy level's fingerprints expected for this species (the
            program's prediction; the target's if it cannot predict)."""
            got = self._expected_fingerprints(
                configurations[can], mc, "energy", version=version
            )
            if got is None:
                got = [fingerprint] if fingerprint else []
            return got

        to_run = []
        for can in order:
            if can in high:
                continue
            fingerprints = expected(can)
            configuration = configurations[can]
            # Already on this configuration in the job (an earlier iteration)?
            E = self._energy_level_from_job(configuration, level, fingerprints)
            if E is not None:
                high[can] = (E, None, "job")
                continue
            row = None
            if store is not None:
                geometry_layer = getattr(self, "_store_layer", {}).get(can, "local")
                kwargs = {"prefer": geometry_layer} if hasattr(store, "lowest") else {}
                row = store.energy(
                    hash_geometry(
                        configuration.atoms.symbols,
                        configuration.atoms.get_coordinates(fractionals=False),
                        configuration.atoms.atomic_numbers,
                    ),
                    energy_fingerprint=fingerprints,
                    model_chemistry=level,
                    **kwargs,
                )
            if row is not None:
                how = self._store_source(store, row)
                if row.get("layer", "local") != geometry_layer and hasattr(
                    store, "layer_name"
                ):
                    # R5: the same geometry, energies from the other layer
                    how += (
                        f"; geometry {store.layer_name(geometry_layer)}, energy "
                        f"{store.layer_name(row['layer'])}"
                    )
                high[can] = (row["E"], None, how)
            else:
                to_run.append(can)
        for can, result in run(to_run).items():
            high[can] = (result.energy, record_of(result), "computed")
            computed_parts[can] = (result.data or {}).get("energy parts")
        # The self-check of the prediction, on what was computed
        for can, (_, record, how) in high.items():
            if how == "computed":
                self._check_prediction(
                    can, configurations[can], "energy", record, mc=mc
                )

        for can in order:
            E, record, how = high[can]
            low = energies[can]
            combined = {"E": E, "E geometry level": low["E"]}
            if "ZPE" in low:
                combined["ZPE"] = low["ZPE"]
                combined["H0"] = E + low["ZPE"]
            if "H298" in low:
                combined["H298"] = E + (low["H298"] - low["E"])
            if "imaginary" in low:
                combined["imaginary"] = low["imaginary"]
            energies[can] = combined
            sources[can] = f"{sources.get(can, '')}; energy {how}".lstrip("; ")
            if how == "computed":
                self._store_energy_level(
                    configurations[can], level, E, record, computed_parts.get(can)
                )
        printer.important(
            __(
                f"\n{level} energies: "
                + ", ".join(f"{can} {high[can][2]}" for can in order)
                + ".",
                indent=4 * " ",
                wrap=True,
                dedent=False,
            )
        )

    def _composite_energies(self, mc, store, order, configurations, energies, sources):
        """The energies of a composite that combines several calculations, on the
        species' geometries: one CompositeEvaluator over all the species.

        Each component is looked up first, never the composite: in the job (an
        earlier iteration's property at the component's level) and in the fragment
        store's energy table by the component's own fingerprint. The rest run, each
        component as one batch of its program. Each computed component is stored
        as its own energy property (so the store's harvest keeps components only),
        and the composite's record, E0 and fingerprint under ``composite ...``
        names, which the harvest does not read."""
        from composite_step.evaluate import CompositeEvaluator
        from seamm_thermochemistry.geometry_store import hash_geometry

        level = mc["level"]
        recipe = getattr(self, "_energy_recipe", None)
        directory = self._species_directory()
        directory = directory / "composite" if directory is not None else None
        keys = {can: f"species_{i}" for i, can in enumerate(order)}
        by_key = {v: k for k, v in keys.items()}

        def lookup(source, component_mc, key, configuration):
            can = by_key[key]
            component = component_mc["level"]
            fingerprints = (
                self._expected_fingerprints(configuration, component_mc, "energy") or []
            )
            E = self._energy_level_from_job(configuration, component, fingerprints)
            if E is not None:
                record = {}
                for key, base in (
                    ("fingerprint", "settings fingerprint"),
                    ("energy parts", "energy parts"),
                ):
                    name = f"{base}#Thermochemistry#{component}"
                    if configuration.properties.exists(name):
                        record[key] = configuration.properties.get(name)[name]["value"]
                return E, "job", record
            if store is None:
                return None
            geometry_layer = getattr(self, "_store_layer", {}).get(can, "local")
            kwargs = {"prefer": geometry_layer} if hasattr(store, "lowest") else {}
            row = store.energy(
                hash_geometry(
                    configuration.atoms.symbols,
                    configuration.atoms.get_coordinates(fractionals=False),
                    configuration.atoms.atomic_numbers,
                ),
                energy_fingerprint=fingerprints,
                model_chemistry=component,
                **kwargs,
            )
            if row is None:
                return None
            record = {
                "fingerprint": row.get("energy_fingerprint"),
                "program": row.get("program"),
                "version": row.get("version"),
                "energy parts": row.get("byproducts") or None,
            }
            return row["E"], self._store_source(store, row), record

        results = {}
        with CompositeEvaluator(
            self,
            mc,
            recipe,
            properties=("energy",),
            directory=directory,
            name="Thermochemistry_composite",
            component_lookup=lookup,
        ) as evaluator:
            for can in order:
                evaluator.submit(configurations[can], key=keys[can])
            for result in evaluator.results():
                results[by_key[result.key]] = result
        for can in order:
            result = results.get(can)
            if result is None or not result.ok or result.energy is None:
                reason = result.reason if result is not None else "no result"
                raise RuntimeError(f"The {level} energy of {can} failed: {reason}")

        for can in order:
            result = results[can]
            composite = result.data["composite"]
            record = result.data["calculation"]
            E = result.energy
            low = energies[can]
            combined = {"E": E, "E geometry level": low["E"]}
            if "ZPE" in low:
                combined["ZPE"] = low["ZPE"]
                combined["H0"] = E + low["ZPE"]
            if "H298" in low:
                combined["H298"] = E + (low["H298"] - low["E"])
            if "imaginary" in low:
                combined["imaginary"] = low["imaginary"]
            energies[can] = combined
            hows = ", ".join(
                f"{source} {c['how']}" for source, c in composite["components"].items()
            )
            sources[can] = f"{sources.get(can, '')}; energy {hows}".lstrip("; ")
            for source, c in composite["components"].items():
                if str(c["how"]).startswith("computed"):
                    self._store_energy_level(
                        configurations[can],
                        c["model chemistry"],
                        c["energy"],
                        record["components"][source]["record"],
                        c.get("parts"),
                    )
            self._store_composite(configurations[can], level, E, record)
        printer.important(
            __(
                f"\n{level} energies: "
                + "; ".join(
                    f"{can} ({sources[can].split('energy ', 1)[-1]})" for can in order
                )
                + ".",
                indent=4 * " ",
                wrap=True,
                dedent=False,
            )
        )

    @staticmethod
    def _store_composite(configuration, level, E, record):
        """A composite's E0 (E_h), fingerprint and record as properties
        ``composite energy|fingerprint|calculation#Thermochemistry#<level>``: with
        the results, never as an energy row of the fragment store (whose harvest
        reads ``total energy`` and ``settings fingerprint`` only)."""
        properties = configuration.properties
        values = [
            ("composite energy", Q_(E, "kJ/mol").m_as("E_h"), "float", "E_h"),
            ("composite fingerprint", record.get("fingerprint"), "str", None),
            ("composite calculation", record, "json", None),
        ]
        for base, value, kind, units in values:
            name = f"{base}#Thermochemistry#{level}"
            if not properties.exists(name):
                properties.add(name, kind, units=units, description=base)
            properties.put(name, value)

    @staticmethod
    def _store_energy_level(configuration, level, E, record, parts=None):
        """The energy-level energy (E_h), its settings fingerprint, the
        calculation record and the parts of the energy (as the program reported
        them) as properties ``...#Thermochemistry#<level>``, for the fragment
        store's harvest."""
        properties = configuration.properties
        values = [
            ("total energy", Q_(E, "kJ/mol").m_as("E_h"), "float", "E_h"),
        ]
        if record:
            if record.get("fingerprint"):
                values.append(
                    ("settings fingerprint", record["fingerprint"], "str", None)
                )
            values.append(("calculation", record, "json", None))
        if parts:
            values.append(("energy parts", parts, "json", None))
        for base, value, kind, units in values:
            name = f"{base}#Thermochemistry#{level}"
            if not properties.exists(name):
                properties.add(name, kind, units=units, description=base)
            properties.put(name, value)

    @staticmethod
    def _target_smiles(configuration, charge, multiplicity):
        """The SMILES of the target, checked against the configuration's charge
        and spin multiplicity, which the CBH reactions and the ATcT lookups
        rely on. A mismatch is an error: molsystem renders a radical centre as
        an anion (``[CH2-]C`` for the ethyl radical), which would silently
        give the reactions of the wrong species."""
        from seamm_thermochemistry import cbh

        smiles = configuration.to_smiles(canonical=True, isomeric=True)
        _, q, radicals, _ = cbh._composition(smiles)
        unpaired = multiplicity - 1
        # Fewer unpaired electrons than the SMILES shows is allowed: the structure's
        # multiplicity wins (SiCl2 and BH are singlets, written as diradicals)
        if q != charge or (unpaired - radicals) % 2 != 0:
            raise RuntimeError(
                f"The SMILES of the structure, {smiles}, implies charge {q} and "
                f"{radicals} unpaired electron(s), but the configuration has charge "
                f"{charge} and spin multiplicity {multiplicity}. Check the charge "
                "and multiplicity of the system. (A radical centre rendered as an "
                "anion, e.g. [CH2-]C for the ethyl radical, is a known limitation "
                "of the SMILES from a structure; build the system from the "
                "radical's SMILES.)"
            )
        return smiles

    def _find_fragment(self, system_db, can):
        """The configuration of the system named CBH:<can>, or None."""
        name = f"CBH:{can}"
        for system in system_db.get_systems():
            if system.name == name:
                return system.configuration
        return None

    def _compute_species(
        self, _P, can, generated, system=None, configuration=None, embed=False
    ):
        """Create (if needed) and compute one species with the subflowchart.
        `embed`: start a given configuration from its lowest MMFF conformer."""
        printer.important(
            __(f"\nSpecies {can}:", indent=4 * " ", wrap=False, dedent=False)
        )
        system_db = self.get_variable("_system_db")
        if configuration is not None and embed:
            self._embed_lowest(
                configuration, generated, int(_P.get("conformer embeddings", 20))
            )
        if configuration is None:
            configuration = self._find_fragment(system_db, can)
            if configuration is None:
                system = system_db.create_system(name=f"CBH:{can}", make_current=True)
                configuration = system.create_configuration(
                    name="initial", make_current=True
                )
                configuration.from_smiles(generated, name=f"CBH:{can}")
                configuration.charge = cbh_charge(generated)
                configuration.spin_multiplicity = cbh_multiplicity(generated)
                if str(_P.get("initial conformer", "lowest")).startswith("lowest"):
                    self._embed_lowest(
                        configuration,
                        generated,
                        int(_P.get("conformer embeddings", 20)),
                    )
            else:
                system = configuration.system
        system_db.system = system
        system.configuration = configuration
        self._working_configuration = configuration

        label = can.replace("/", "_").replace("\\", "_")
        self._subdirectory = f"species_{label}"
        try:
            if getattr(self, "_engine_mode", "subflowchart") == "Model Chemistry":
                result = self._compute_species_mc(_P, can, configuration)
            elif _P["thermal corrections"] != "harmonic":
                atoms = self.ase_atoms(configuration)
                E = Q_(atoms.get_potential_energy(), "eV").m_as("kJ/mol")
                printer.important(
                    __(
                        f"    enthalpy of formation (from the engine) = {E:.3f} kJ/mol",
                        indent=4 * " ",
                        wrap=False,
                        dedent=False,
                    )
                )
                result = {"H298": E, "E": E}
            else:
                if _yes(_P["optimize fragments"]) and configuration.n_atoms > 1:
                    self._optimize(configuration, _P)
                self._XYZ0 = np.array(
                    configuration.atoms.get_coordinates(fractionals=False)
                )
                self._E0 = None
                XYZ_save = configuration.atoms.coordinates
                self._vibrations = self.run_ase_Hessian(
                    step_size=_P["step size"].m_as("Å"),
                    on_error=_P["on error"],
                    on_success=_P["on success"],
                )
                configuration.atoms.coordinates = XYZ_save
                data = self.harmonic_thermochemistry(
                    _P, [T_REF], "K", [1.0], "bar", print_frequencies=False
                )
                if "total_energy" not in data:
                    raise RuntimeError(
                        f"No energy for the undisplaced structure of {can}"
                    )
                E = Q_(data["total_energy"], "E_h").m_as("kJ/mol")
                result = {
                    "E": E,
                    "ZPE": data["ZPE"],
                    "H0": E + data["ZPE"],
                    "H298": E + data["H thermal"][0][0],
                }
                result["imaginary"] = data.get("N saddle modes", 0)
                # Keep them on the configuration for "existing results" runs
                self._store_species_properties(configuration, data)
        finally:
            self._subdirectory = None
        try:
            self._write_species_sdf(configuration, can, result)
        except Exception as e:
            # Diagnostic output: it must never fail the species
            logger.warning(f"Could not write the structure of {can}: {e}")
        return result

    def _species_directory(self):
        """The step's directory, or None when it has none (outside a flowchart)."""
        try:
            return Path(self.directory)
        except Exception:
            return None

    def _write_species_sdf(self, configuration, can, result):
        """Write the species' structure, as computed, to ``<SMILES>.sdf`` in the
        step's directory and append it to ``species.sdf``, with its energies as
        SDF data fields, so the conformations can be checked; and print the
        heavy-atom dihedrals of the final structure next to the starting ones."""
        directory = self._species_directory()
        if directory is None:
            return
        label = can.replace("/", "_").replace("\\", "-")
        symbols = configuration.atoms.symbols
        xyz = configuration.atoms.get_coordinates(fractionals=False)
        dihedrals = self._chain_dihedrals(symbols, xyz)
        try:
            text = configuration.to_sdf_text()
        except Exception as e:  # pragma: no cover - needs Open Babel
            logger.warning(f"Could not write the structure of {can}: {e}")
            return
        name = configuration.system.name if configuration.system else ""
        fields = {
            "SMILES": can,
            "role": "fragment" if str(name).startswith("CBH:") else "target",
            "model chemistry": self.model or "",
        }
        for key, title in (
            ("E", "E (kJ/mol)"),
            ("ZPE", "ZPE (kJ/mol)"),
            ("H0", "H(0 K) (kJ/mol)"),
            ("H298", "H(298.15 K) (kJ/mol)"),
        ):
            value = result.get(key)
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                fields[title] = f"{value:.3f}"
        if "imaginary" in result:
            fields["imaginary modes"] = str(result["imaginary"])
        body = text.rstrip()
        if body.endswith("$$$$"):
            body = body[: -len("$$$$")].rstrip()
        # A data item ends with a blank line. Stripping removed the one after the
        # last of SEAMM's own items (its stored properties), so put it back, or a
        # reader joins that value with the next field.
        body += "\n" if body.endswith("M  END") else "\n\n"
        body += "".join(f"> <{k}>\n{v}\n\n" for k, v in fields.items())
        body += "$$$$\n"
        try:
            directory.mkdir(parents=True, exist_ok=True)
            (directory / f"{label}.sdf").write_text(body)
            with open(directory / "species.sdf", "a") as fd:
                fd.write(body)
        except OSError as e:  # pragma: no cover
            logger.warning(f"Could not write the structure of {can}: {e}")
            return
        printer.important(
            __(
                f"    structure: {label}.sdf"
                + (f"; final{dihedrals[1:]}" if dihedrals else ""),
                indent=4 * " ",
                wrap=False,
                dedent=False,
            )
        )

    def _compute_species_mc(self, _P, can, configuration):
        """One species through the Model Chemistry: energy only, or optimize
        + Hessian + harmonic thermochemistry."""
        if _P["thermal corrections"] != "harmonic":
            E = self.optimize_and_energy(
                configuration,
                optimize=(_yes(_P["optimize fragments"])),
                fmax=_P["fragment convergence"].m_as("eV/Å"),
                max_steps=int(_P["fragment max steps"]),
            )
            printer.important(
                __(
                    f"    enthalpy of formation (from the engine) = {E:.3f} kJ/mol",
                    indent=4 * " ",
                    wrap=False,
                    dedent=False,
                )
            )
            self._store_fingerprints(configuration, None, can)
            return {"H298": E, "E": E}
        delta_bohr = Q_(_P["step size"].m_as("Å"), "Å").m_as("bohr")
        E0, hessian, how = self.optimize_and_hessian(
            configuration,
            delta_bohr,
            optimize=(_yes(_P["optimize fragments"])),
            fmax=_P["fragment convergence"].m_as("eV/Å"),
            max_steps=int(_P["fragment max steps"]),
        )
        from ase.vibrations import VibrationsData

        atoms = self.ase_atoms(configuration)
        n = len(atoms)
        hessian = Q_(np.asarray(hessian), "hartree/bohr**2").m_as("eV/Å**2")
        self._vibrations = VibrationsData(atoms, hessian.reshape(n, 3, n, 3))
        self._E0 = Q_(E0, "kJ/mol").m_as("eV")
        printer.important(
            __(f"    Hessian: {how}", indent=4 * " ", wrap=False, dedent=False)
        )
        data = self.harmonic_thermochemistry(
            _P, [T_REF], "K", [1.0], "bar", print_frequencies=False
        )
        E = Q_(data["total_energy"], "E_h").m_as("kJ/mol")
        result = {
            "E": E,
            "ZPE": data["ZPE"],
            "H0": E + data["ZPE"],
            "H298": E + data["H thermal"][0][0],
            "imaginary": data.get("N saddle modes", 0),
        }
        self._store_species_properties(configuration, data)
        self._store_fingerprints(configuration, how, can)
        return result

    # -- fingerprints of each species (P5a of the fragment store) ---------------
    #
    # A program's settings can depend on the species, not only on the model
    # chemistry: ORCA uses exact exchange (NoCOSX) for lone atoms and for Na, Mg,
    # Zn, B and P (orca_step#44), so a boron species' revDSD energy has another
    # fingerprint than methane's. Each species is therefore looked up -- in the
    # job's database and in the fragment store -- with the fingerprints the
    # program predicts for IT (seamm_exec.expected_fingerprints). The target's
    # fingerprints are the fallback only for programs that cannot predict.

    def _reset_fingerprints(self):
        """Forget the previous run's fingerprints (a Loop runs this node again)."""
        self._settings_fingerprint = None
        self._frequency_fingerprint = None
        self._fingerprint_version = None
        self._hessian_how = None

    def _geometry_mc(self):
        """The geometry level's ``_model_chemistry`` wrapper, or None (e.g. the
        subflowchart route, or no flowchart variables)."""
        try:
            return self.model_chemistry()
        except Exception:
            return None

    def _expected_fingerprints(self, structure, mc, kind, version=None):
        """The program's predicted fingerprints of a structure at a model
        chemistry ([str]), or None when it cannot predict them."""
        if mc is None or structure is None:
            return None
        # The program's version is hashed: without one from a calculation of
        # this run (e.g. "existing results" computes nothing) a prediction
        # would be wrong, so make none
        version = version or getattr(self, "_fingerprint_version", None)
        if not version:
            return None
        try:
            from seamm_exec import expected_fingerprints
        except ImportError:  # pragma: no cover - an older seamm-exec
            return None
        try:
            provider = self.flowchart.plugin_manager.get(mc["step"])
        except Exception:
            return None
        if provider is None:
            return None
        return expected_fingerprints(
            provider,
            mc,
            structure,
            kind=kind,
            properties=("energy",),
            version=version,
        )

    @staticmethod
    def _species_structure(can):
        """A stand-in structure for a species not yet built: its elements (with
        the hydrogens), charge and multiplicity, which are all a program's
        settings depend on; the coordinates are not settings."""
        from rdkit import Chem
        from seamm_exec import Geometry

        mol = Chem.MolFromSmiles(can)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
        return Geometry(
            numbers,
            [[0.0, 0.0, 0.0]] * len(numbers),
            charge=cbh_charge(can),
            multiplicity=cbh_multiplicity(can),
        )

    def _geometry_level_fingerprints(self, can, configuration=None):
        """(energy fingerprints, frequency fingerprints) expected for a species
        at the geometry level, each a list or None. The frequency ones are the
        analytic Hessian's when the target's Hessian was analytic, else the
        energy's (finite differences of gradients)."""
        mc = self._geometry_mc()
        structure = configuration
        if structure is None:
            try:
                structure = self._species_structure(can)
            except Exception:  # pragma: no cover - RDKit missing or odd SMILES
                structure = None
        energy = self._expected_fingerprints(structure, mc, "energy")
        if getattr(self, "_hessian_how", None) == "analytic":
            frequency = self._expected_fingerprints(structure, mc, "hessian")
        else:
            frequency = energy
        if energy is None and getattr(self, "_settings_fingerprint", None):
            energy = [self._settings_fingerprint]
        if frequency is None and getattr(self, "_frequency_fingerprint", None):
            frequency = [self._frequency_fingerprint]
        return energy, frequency

    def _check_prediction(self, can, configuration, kind, record, mc=None):
        """The self-check: warn when a computed species' fingerprint is not one
        the program predicted, i.e. the prediction has drifted from the
        program's rules. The species is still used and stored with its own
        fingerprint; a drift costs reuse, never correctness."""
        if not record or not record.get("fingerprint"):
            return
        if mc is None:
            mc = self._geometry_mc()
        expected = self._expected_fingerprints(
            configuration, mc, kind, version=record.get("version")
        )
        if expected is None or record["fingerprint"] in expected:
            return
        text = (
            f"    Warning: the {kind} fingerprint of {can} "
            f"({record['fingerprint'][:15]}...) is not the one predicted for it "
            f"({', '.join(f[:15] + '...' for f in expected)}), so it will not be "
            "found in the fragment store by later jobs. Please report this."
        )
        logger.warning(text.strip())
        printer.important(__(text, indent=4 * " ", wrap=True, dedent=False))

    @staticmethod
    def _property(configuration, name, units=None):
        """A property's value (in `units` if given), or None."""
        properties = configuration.properties
        try:
            if not properties.exists(name):
                return None
            values = properties.get(pattern=name, match="exact")
        except Exception:
            return None
        if not values or name not in values:
            return None
        value = values[name]
        if isinstance(value, dict) and "value" in value:
            value = value["value"]
        if units is not None and value is not None:
            stored = properties.units(name)
            value = Q_(float(value), stored if stored else units).m_as(units)
        return value

    def _from_job(self, _P, system_db, can, split):
        """A species already computed in this job (a system named CBH:<can>,
        e.g. by an earlier iteration of a Loop), with the fingerprints expected
        for it at the geometry level, as ({"E", "ZPE", "H0", "H298"},
        configuration), or None.

        `split`: whether the energies are at another level, so that the
        frequency fingerprint is the key; otherwise the settings fingerprint.
        A configuration with an imaginary mode above 20i is not reused.
        """
        from seamm_thermochemistry.geometry_store import IMAGINARY_TOLERANCE

        energy_fps, frequency_fps = self._geometry_level_fingerprints(can)
        expected = frequency_fps if split else energy_fps
        if not expected:
            return None
        harmonic = _P["thermal corrections"] == "harmonic"
        model = self.model or "unknown"
        key = "frequency fingerprint" if split else "settings fingerprint"
        for system in system_db.get_systems():
            if system.name != f"CBH:{can}":
                continue
            current = system.configuration
            candidates = [current] + [
                c for c in system.configurations if c.id != current.id
            ]
            for configuration in candidates:
                fingerprint = self._property(
                    configuration, f"{key}#Thermochemistry#{model}"
                )
                if fingerprint not in expected:
                    continue
                E = self._property(
                    configuration, f"total energy#Thermochemistry#{model}", "kJ/mol"
                )
                if E is None:
                    continue
                result = {"E": E}
                if harmonic:
                    ZPE = self._property(
                        configuration, f"ZPE#Thermochemistry#{model}", "kJ/mol"
                    )
                    H = self._property(
                        configuration, f"H(298.15 K)#Thermochemistry#{model}", "kJ/mol"
                    )
                    if ZPE is None or H is None:
                        continue
                    frequencies = self._property(
                        configuration, f"frequencies#Thermochemistry#{model}"
                    )
                    if frequencies and any(
                        -float(f) > IMAGINARY_TOLERANCE for f in frequencies
                    ):
                        continue
                    result.update(ZPE=ZPE, H0=E + ZPE, H298=H)
                else:
                    result["H298"] = E
                return result, configuration
        return None

    def _energy_level_from_job(self, configuration, level, expected):
        """The energy level's energy (kJ/mol) already on this configuration with
        one of the expected fingerprints, or None.

        The calculation must have run on the configuration's current geometry
        (its record's geometry fingerprint): a species computed again in place
        has new coordinates but keeps the old energy-level properties."""
        from seamm_thermochemistry.geometry_store import hash_geometry

        if not expected:
            return None
        fingerprint = self._property(
            configuration, f"settings fingerprint#Thermochemistry#{level}"
        )
        if fingerprint not in expected:
            return None
        record = self._property(configuration, f"calculation#Thermochemistry#{level}")
        if isinstance(record, str):
            try:
                record = json.loads(record)
            except ValueError:
                record = None
        if not isinstance(record, dict) or not record.get("geometry fingerprint"):
            return None
        here = hash_geometry(
            configuration.atoms.symbols,
            configuration.atoms.get_coordinates(fractionals=False),
            configuration.atoms.atomic_numbers,
        )
        if record["geometry fingerprint"] != here:
            return None
        return self._property(
            configuration, f"total energy#Thermochemistry#{level}", "kJ/mol"
        )

    def _store_fingerprints(self, configuration, how, can=None):
        """Store the engine's settings fingerprints as the properties
        ``settings fingerprint``, ``element fingerprints`` and ``calculation``
        (and the ``frequency ...`` three for a Hessian: the analytic-Hessian
        record, or the energy record for finite differences) under
        ``#Thermochemistry#<model>``, the names the fragment store's harvest
        reads. Nothing is stored when the engine gives no fingerprints.

        The first species of a run (the target) also sets the fallback
        fingerprints, the program's version and the Hessian route; each
        species' fingerprints are checked against the predicted ones."""
        records = getattr(self, "_fingerprints", None)
        if not records:
            return
        energy = records.get("energy")
        if energy and energy.get("fingerprint"):
            if not getattr(self, "_settings_fingerprint", None):
                self._settings_fingerprint = energy["fingerprint"]
            if not getattr(self, "_fingerprint_version", None):
                self._fingerprint_version = energy.get("version")
        if how is not None and not getattr(self, "_hessian_how", None):
            self._hessian_how = how
        label = can or "the species"
        self._check_prediction(label, configuration, "energy", energy)
        kinds = [("", energy)]
        if how is not None:
            frequency = records.get("hessian") if how == "analytic" else energy
            if frequency and frequency.get("fingerprint"):
                if not getattr(self, "_frequency_fingerprint", None):
                    self._frequency_fingerprint = frequency["fingerprint"]
                if how == "analytic":
                    self._check_prediction(label, configuration, "hessian", frequency)
            kinds.append(
                ("frequency ", records.get("hessian") if how == "analytic" else energy)
            )
        model = self.model or "unknown"
        properties = configuration.properties
        for prefix, record in kinds:
            if not record or not record.get("fingerprint"):
                continue
            for name, value, kind in (
                (f"{prefix}settings fingerprint", record["fingerprint"], "str"),
                (
                    f"{prefix}element fingerprints",
                    record.get("element fingerprints", {}),
                    "json",
                ),
                (f"{prefix}calculation", record, "json"),
            ):
                if prefix and name.startswith("frequency settings"):
                    name = "frequency fingerprint"
                full = f"{name}#Thermochemistry#{model}"
                if not properties.exists(full):
                    properties.add(full, kind, units=None, description=name)
                properties.put(full, value)

    def _optimize(self, configuration, _P):
        """Optimize the structure with ASE's BFGS through the subflowchart."""
        from ase.optimize import BFGS

        atoms = self.ase_atoms(configuration)
        fmax = _P["fragment convergence"].m_as("eV/Å")
        optimizer = BFGS(atoms, logfile=None)
        converged = optimizer.run(fmax=fmax, steps=int(_P["fragment max steps"]))
        configuration.atoms.set_coordinates(atoms.positions, fractionals=False)
        if not converged:
            printer.important(
                __(
                    f"    Warning: the optimization did not converge in "
                    f"{_P['fragment max steps']} steps.",
                    indent=4 * " ",
                    wrap=False,
                    dedent=False,
                )
            )

    def _store_species_properties(self, configuration, data):
        """Store E, ZPE and H(298.15 K) on the configuration."""
        model = self.model or "unknown"
        properties = configuration.properties
        for name, value, units, description in (
            (f"total energy#Thermochemistry#{model}", data["total_energy"], "E_h", "E"),
            (f"ZPE#Thermochemistry#{model}", data["ZPE"], "kJ/mol", "ZPE"),
            (
                f"H(0 K)#Thermochemistry#{model}",
                Q_(data["total_energy"], "E_h").m_as("kJ/mol") + data["ZPE"],
                "kJ/mol",
                "H(0 K)",
            ),
            (
                f"H(298.15 K)#Thermochemistry#{model}",
                Q_(data["total_energy"], "E_h").m_as("kJ/mol")
                + data["H thermal"][0][0],
                "kJ/mol",
                "H(298.15 K)",
            ),
        ):
            if not properties.exists(name):
                properties.add(name, "float", units=units, description=description)
            properties.put(name, value)
        # For the fragment store's geometry rows: the frequencies (imaginary ones
        # negative), the rotational symmetry number and the moments of inertia
        extras = [
            ("frequencies", data.get("vibrational frequencies"), "json", "cm^-1"),
            ("symmetry number", data.get("symmetry number"), "int", None),
        ]
        try:
            moments = self.ase_atoms(configuration).get_moments_of_inertia().tolist()
            extras.append(("moments of inertia", moments, "json", "amu*Å^2"))
        except Exception:  # pragma: no cover
            pass
        for base, value, kind, units in extras:
            if value is None:
                continue
            name = f"{base}#Thermochemistry#{model}"
            if kind == "json":
                value = [float(v) for v in value]
            if not properties.exists(name):
                properties.add(name, kind, units=units, description=base)
            properties.put(name, value)

    def _energies_from_properties(self, _P, configuration, can):
        """E, ZPE and H(298.15 K) from the configuration's properties."""
        properties = configuration.properties
        result = {}

        def get(name):
            if name in ("", "none"):
                return None
            values = properties.get(pattern=name, match="glob")
            if not values:
                return None
            if len(values) > 1:
                raise RuntimeError(
                    f"Several properties of {can} match '{name}': {sorted(values)}"
                )
            (key, value), *_ = values.items()
            if isinstance(value, dict):  # {"value": ..., "sid": ..., "cid": ...}
                value = value["value"]
            units = properties.units(key)
            return Q_(float(value), units if units else "kJ/mol").m_as("kJ/mol")

        E = get(_P["energy property"])
        if E is None:
            raise RuntimeError(
                f"No property matching '{_P['energy property']}' on the system for "
                f"{can}."
            )
        result["E"] = E
        if _P["thermal corrections"] != "harmonic":
            result["H298"] = E
            return result
        ZPE = get(_P["ZPE property"])
        H = get(_P["enthalpy property"])
        E_low = get(_P["low-level energy property"])
        if E_low is None:
            E_low = E
        if ZPE is not None:
            result["ZPE"] = ZPE
            result["H0"] = E + ZPE
        if H is not None:
            result["H298"] = E + (H - E_low)
        return result

    def _cite_cbh(self, tables):
        for key, alias, note in (
            ("Ruscic2004", "ATcT", "The Active Thermochemical Tables."),
            ("Ruscic2005", "ATcT2005", "The Active Thermochemical Tables."),
            ("Ramabhadran2011", "CBH", "The connectivity-based hierarchy."),
            (
                "Raghavachari2023",
                "CBH2023",
                "The graph formulation of the connectivity-based hierarchy.",
            ),
        ):
            if key in self._bibliography:
                self.references.cite(
                    raw=self._bibliography[key],
                    alias=alias,
                    module="thermochemistry_step",
                    level=1,
                    note=note,
                )
        version = tables.version
        raw = (
            "@misc{ATcTvalues,\n"
            "    author = {Ruscic, Branko and Bross, David H.},\n"
            "    title = {Active Thermochemical Tables (ATcT) values based on ver. "
            f"{version} of the Thermochemical Network}},\n"
            "    publisher = {Argonne National Laboratory},\n"
            "    address = {Lemont, Illinois, USA},\n"
            "    howpublished = {available at ATcT.anl.gov},\n"
            "}"
        )
        self.references.cite(
            raw=raw,
            alias=f"ATcT-{version}",
            module="thermochemistry_step",
            level=1,
            note="The ATcT values used as references.",
        )


def rung_spread(results, usable):
    """(spread, rungs): the range of ΔfH°(298.15 K) over the usable rungs from
    CBH-1 up, and those rungs; (None, rungs) with fewer than two."""
    rungs = [r for r in usable if r >= 1 and "DfH 298.15 K" in results[r]]
    if len(rungs) < 2:
        return None, rungs
    values = [results[r]["DfH 298.15 K"] for r in rungs]
    return max(values) - min(values), rungs


def species_only(P):
    """Whether the CBH approach computes only the target, for the store."""
    return str(P.get("cbh calculation", "enthalpy of formation")).startswith(
        "the species only"
    )


def _yes(value):
    """A yes/no parameter, which arrives as a bool or a string."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("yes", "true", "1", "on")


def cbh_charge(smiles):
    from seamm_thermochemistry.cbh import _composition

    return _composition(smiles)[1]


def cbh_multiplicity(smiles):
    """The spin multiplicity of a CBH species: its ground state, from the tables of
    singlets and triplets that a SMILES does not show (SiCl2 and BH are singlets,
    O2 and OCCO triplets), else its unpaired electrons + 1."""
    from seamm_thermochemistry import cbh

    return cbh.ground_state_multiplicity(smiles)


def _unused():  # pragma: no cover
    return Counter, math


def _reference_cell(entry, conversion, temperature, column):
    """A reference value for the species table: blank if none, starred if
    converted from the other temperature with the computed H298 - H0."""
    if not entry:
        return ""
    if conversion is not None and conversion[0] == temperature:
        return f"{conversion[1]:.3f}*"
    value = entry.get(column)
    return "" if value is None else f"{value:.3f}"
