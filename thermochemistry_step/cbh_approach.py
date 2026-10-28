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
import logging
import math
import os

import numpy as np
from tabulate import tabulate

from seamm_util import Q_, parse_list
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
printer = printing.getPrinter("Thermochemistry")

T_REF = 298.15

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

        system_db = self.get_variable("_system_db")
        target_system, target_configuration = self.get_system_configuration()
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

        text = f"The target is {smiles}"
        if charge != 0 or multiplicity != 1:
            text += f" with charge {charge} and multiplicity {multiplicity}"
        text += ".\n"

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
        tables = ATcT(version=version, base=base)
        text += (
            f"Reference values from ATcT ver. {tables.version} of the "
            f"Thermochemical Network ({tables.directory}).\n"
        )

        def lookup(can):
            _, q, radicals, _ = cbh._composition(can)
            return tables.lookup(can, charge=q, multiplicity=radicals + 1)

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

        def reference(can):
            """The reference entry under the policy: an ATcT entry (with
            atct_coefficients {id: 1}) or the ladder's entry."""
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
        if _yes(_P["compare with ATcT"]):
            try:
                target_entry = tables.lookup(
                    target, charge=charge, multiplicity=multiplicity
                )
            except (MissingEntry, AmbiguousEntry):
                target_entry = None

        # Report the reactions
        lines = []
        for rung in rungs:
            if rung in unavailable:
                lines.append(f"CBH-{rung}: not available. {unavailable[rung]}")
            elif reactions[rung].trivial:
                lines.append(f"CBH-{rung}: trivial, the rung reproduces the target.")
            elif rung in missing:
                lines.append(
                    f"{reactions[rung]}\n    not usable: no ATcT entry for "
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
        store = self._fragment_store(_P)
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
                _P, target, smiles, target_system, target_configuration
            )
            sources[target] = "computed"
            for can, generated in fragments.items():
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
            # Back to the target
            system_db.system = target_system
            target_system.configuration = target_configuration
            self._working_configuration = target_configuration

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

        results = {}
        for rung in usable:
            reaction = reactions[rung]
            result = {"rung": rung, "reaction": str(reaction)}
            for key, label in (("H0", "0 K"), ("H298", "298.15 K")):
                if all(key in energies[s] for s in reaction.coefficients):
                    drH = sum(
                        c * energies[s][key] for s, c in reaction.coefficients.items()
                    )
                    column = "dfH0_0K" if key == "H0" else "dfH0_298K"
                    ref = sum(
                        c * entries[s][column] for s, c in reaction.references.items()
                    )
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
            covariance = (
                tables.covariance(ids, fetch=(_P["covariances"] == "fetch if needed"))
                if ids
                else None
            )
            u_atct, used = reaction_uncertainty(coefficients, uncertainties, covariance)
            u_method = 0.0
            if method_c2 > 0 and ladder is not None:
                rms = ladder.method_uncertainty.get("rms_0K")
                if rms:
                    u_method = rms * method_c2**0.5
            u = (u_atct**2 + u_method**2) ** 0.5
            result["uncertainty"] = u
            result["uncertainty ATcT"] = u_atct
            result["uncertainty method"] = (
                "ATcT covariance" if used else "quadrature (no covariance cached)"
            )
            if u_method:
                result["uncertainty ladder method"] = u_method
                result["uncertainty method"] += " + ladder method term"
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
                    f"{e['ZPE']:.3f}" if "ZPE" in e else "",
                    f"{e['H298'] - e['E']:.3f}" if "E" in e and "H298" in e else "",
                    (
                        (
                            f"ladder (CBH-{entry['rung']})"
                            if entry.get("source") == "laddered"
                            else entry["id"]
                        )
                        if entry
                        else ""
                    ),
                    (
                        f"{entry['dfH0_0K']:.3f}"
                        if entry and entry["dfH0_0K"] is not None
                        else ""
                    ),
                    f"{entry['dfH0_298K']:.3f}" if entry else "",
                    (
                        "exact"
                        if entry and entry["exact"]
                        else f"{entry['uncertainty_298K']:.3f}" if entry else ""
                    ),
                ]
            )
        text = "\nThe species (energies in kJ/mol):\n"
        text += tabulate(
            species_rows,
            headers=[
                "Species",
                "Source",
                "E",
                "ZPE",
                "H(298.15 K) - E",
                "Reference",
                "ATcT ΔfH°(0 K)",
                "ATcT ΔfH°(298.15 K)",
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
        if any(
            r["uncertainty method"].startswith("quadrature") for r in results.values()
        ):
            text += (
                "\nNo ATcT covariance block was cached for some reactions, so their "
                "reference uncertainties ignore the correlations (quadrature); run "
                "'seamm-thermochemistry-atct covariance' to cache them."
            )
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        # The results: the highest usable rung
        best = results[usable[-1]]
        data["CBH rung"] = best["rung"]
        data["CBH reaction"] = best["reaction"]
        if "DfH 0 K" in best:
            data["DfH 0K"] = best["DfH 0 K"]
        if "DfH 298.15 K" in best:
            data["DfH 298K"] = best["DfH 298.15 K"]
        data["DfH uncertainty"] = best["uncertainty"]
        data["ATcT version"] = tables.version
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

    # -----------------------------------------------------------------------

    def _fragment_store(self, _P):
        """The fragment store to read, or None if ignored or unreadable."""
        if _P.get("fragment store", "ignore") != "use if available":
            return None
        from seamm_thermochemistry import FragmentStore

        path = _P.get("fragment store file", "default")
        path = None if path in ("default", "") else os.path.expanduser(path)
        try:
            return FragmentStore(path)
        except Exception as e:  # pragma: no cover
            logger.warning(f"The fragment store could not be opened: {e}")
            return None

    def _energies_from_store(self, _P, store, can):
        """The energies of a species from the store at this model chemistry,
        as ({"E", "ZPE", "H0", "H298"}, source), or None. A record must have
        what the thermal-correction mode needs: E and H(298.15 K) with
        harmonic corrections (ZPE too for the 0 K result), E alone without."""
        if store is None:
            return None
        from thermochemistry_step.cbh_approach import cbh_charge, cbh_multiplicity

        model = self.model or "unknown"
        record = store.get(
            can,
            charge=cbh_charge(can),
            multiplicity=cbh_multiplicity(can),
            energy_fingerprint=getattr(self, "_settings_fingerprint", None),
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
        source = "store" if record.get("verified") else "store (unverified)"
        return result, source

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
        if q != charge or unpaired < radicals or (unpaired - radicals) % 2 != 0:
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

    def _compute_species(self, _P, can, generated, system=None, configuration=None):
        """Create (if needed) and compute one species with the subflowchart."""
        system_db = self.get_variable("_system_db")
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
            else:
                system = configuration.system
        system_db.system = system
        system.configuration = configuration
        self._working_configuration = configuration

        label = can.replace("/", "_").replace("\\", "_")
        self._subdirectory = f"species_{label}"
        printer.important(
            __(f"\nSpecies {can}:", indent=4 * " ", wrap=False, dedent=False)
        )
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
                # Keep them on the configuration for "existing results" runs
                self._store_species_properties(configuration, data)
        finally:
            self._subdirectory = None
        return result

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
            self._store_fingerprints(configuration, None)
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
        }
        self._store_species_properties(configuration, data)
        self._store_fingerprints(configuration, how)
        return result

    def _store_fingerprints(self, configuration, how):
        """Store the engine's settings fingerprints as the properties
        ``settings fingerprint``, ``element fingerprints`` and ``calculation``
        (and the ``frequency ...`` three for a Hessian: the analytic-Hessian
        record, or the energy record for finite differences) under
        ``#Thermochemistry#<model>``, the names the fragment store's harvest
        reads. Nothing is stored when the engine gives no fingerprints."""
        records = getattr(self, "_fingerprints", None)
        if not records:
            return
        energy = records.get("energy")
        if energy and energy.get("fingerprint"):
            self._settings_fingerprint = energy["fingerprint"]
        kinds = [("", energy)]
        if how is not None:
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


def _yes(value):
    """A yes/no parameter, which arrives as a bool or a string."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("yes", "true", "1", "on")


def cbh_charge(smiles):
    from seamm_thermochemistry.cbh import _composition

    return _composition(smiles)[1]


def cbh_multiplicity(smiles):
    from seamm_thermochemistry.cbh import _composition

    return _composition(smiles)[2] + 1


def _unused():  # pragma: no cover
    return Counter, math
