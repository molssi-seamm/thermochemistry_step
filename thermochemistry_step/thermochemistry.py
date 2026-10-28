# -*- coding: utf-8 -*-

"""Non-graphical part of the Thermochemistry step in a SEAMM flowchart"""

import importlib
import logging
import pprint  # noqa: F401
import sys
import time

from ase.thermochemistry import IdealGasThermo
from ase.vibrations import VibrationsData
import numpy as np
from tabulate import tabulate

import thermochemistry_step
from thermochemistry_step.hessian import ModelChemistryHessian
from thermochemistry_step.symmetry import rotational_symmetry_number
import molsystem
import seamm
from seamm_ase import ASE_mixin
from seamm_util import getParser, parse_list, Q_, units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Thermochemistry")

# Add this module's properties to the standard properties
path = importlib.resources.files("thermochemistry_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class Thermochemistry(seamm.Node, ASE_mixin, ModelChemistryHessian):
    """
    The non-graphical part of a Thermochemistry step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : ThermochemistryParameters
        The control parameters for Thermochemistry.

    See Also
    --------
    TkThermochemistry,
    Thermochemistry, ThermochemistryParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Thermochemistry",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for thermochemistry in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Thermochemistry {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Thermochemistry", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Thermochemistry",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = thermochemistry_step.metadata
        self.parameters = thermochemistry_step.ThermochemistryParameters()
        self._results = {}
        self._vibrations = None
        self._E0 = None  # energy of the undisplaced structure, eV
        self._XYZ0 = None  # its coordinates

        # Variable used by seamm-ase
        self._data = {}
        self._file_handler = None
        self._last_coordinates = None
        self._logfile = None
        self._step = None

    @property
    def version(self):
        """The semantic version of this module."""
        return thermochemistry_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return thermochemistry_step.__git_revision__

    def model_chemistry_vibrations(self, configuration, step_size):
        """The energy (eV) and an ASE VibrationsData from the Model Chemistry.

        Parameters
        ----------
        configuration : molsystem._Configuration
        step_size : float
            The finite-difference step, Å, if the Hessian is not analytic.

        Returns
        -------
        (float, ase.vibrations.VibrationsData, str)
        """
        delta_bohr = Q_(step_size, "Å").m_as("bohr")
        E0, hessian, how = self.energy_and_hessian(configuration, delta_bohr)
        atoms = self.ase_atoms(configuration)
        hessian = Q_(np.asarray(hessian), "hartree/bohr**2").m_as("eV/Å**2")
        n = len(atoms)
        vibrations = VibrationsData(atoms, hessian.reshape(n, 3, n, 3))
        return Q_(E0, "kJ/mol").m_as("eV"), vibrations, how

    def _subflowchart_model(self, flowchart=None):
        """The model chemistry of the engine in the subflowchart, for naming the
        properties: the last model set by any step, searching sub-subflowcharts."""
        if flowchart is None:
            flowchart = self.subflowchart
        model = None
        node = flowchart.get_node("1").next()
        while node is not None:
            if getattr(node, "model", None) is not None:
                model = node.model
            if getattr(node, "subflowchart", None) is not None:
                tmp = self._subflowchart_model(node.subflowchart)
                if tmp is not None:
                    model = tmp
            node = node.next()
        return model

    @staticmethod
    def _values_and_units(value, default_units):
        """The list of values and the units of a T or P parameter, which may be
        a (string, units) tuple, a quantity, a plain string or a number."""
        if isinstance(value, tuple):
            values, units = value
            return parse_list(values), units
        if hasattr(value, "magnitude"):
            return [value.magnitude], str(value.u)
        if isinstance(value, str):
            return parse_list(value), default_units
        return [float(value)], default_units

    def ase_calculator(self, calculator, *args, **kwargs):
        """Run the subflowchart for one geometry, keeping the energy of the
        undisplaced structure for the absolute thermodynamic functions."""
        result = super().ase_calculator(calculator, *args, **kwargs)
        if self._E0 is None and self._XYZ0 is not None:
            if np.allclose(calculator.atoms.positions, self._XYZ0, atol=1.0e-8):
                # eV, as the calculator reports it
                self._E0 = calculator.results["energy"]
        return result

    def analyze(self, _P, indent="", **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        text = ""
        atoms = self.ase_atoms(self._working_configuration)

        # Get the moments of inertia to see if linear, etc.
        n_atoms = self._working_configuration.n_atoms

        if n_atoms == 1:
            geometry = "monatomic"
            n_zeroes = 3
        else:
            moments = atoms.get_moments_of_inertia().tolist()
            if moments[0] < 0.01 and abs(moments[2] - moments[1]) < 0.01:
                geometry = "linear"
                n_zeroes = 5
            else:
                geometry = "nonlinear"
                n_zeroes = 6

        n_imaginary = 0
        imaginary_frequencies = []
        real_frequencies = []
        real_energies = np.array([])
        if n_atoms > 1:
            # Ignore any imaginary modes in the thermochemistry
            energies = self._vibrations.get_energies()

            n_imaginary = _P["# imaginary modes"]

            # How many imaginary frequencies?
            if isinstance(n_imaginary, int):
                pass
            elif n_imaginary == "transition state":
                n_imaginary = 1
            elif n_imaginary == "perceive":
                vib_text = "\n".join(
                    [
                        f"\t{r:9.4f}, {i:9.4f}i"
                        for r, i in zip(energies.real, energies.imag)
                    ]
                )

                magnitudes = np.abs(energies)
                smallest = sorted(magnitudes)[n_zeroes]
                _type = "imaginary"
                n_imaginary = 0
                for magnitude, is_complex in zip(magnitudes, np.iscomplex(energies)):
                    if _type == "imaginary":
                        if magnitude >= smallest:
                            if not is_complex:
                                raise ValueError(
                                    f"Nominally imaginary mode energy {magnitude:.5f} "
                                    f"is not imaginary. {smallest=:.4f}. This "
                                    "shouldn't happen! Please report this."
                                    "\n\nEnergies:\n" + vib_text
                                )
                            n_imaginary += 1
                        else:
                            _type = "zero"
                            _n_zeroes = 1
                    elif _type == "zero":
                        if magnitude >= smallest:
                            _type = "real"
                        else:
                            _n_zeroes += 1
                if _n_zeroes != n_zeroes:
                    raise ValueError(
                        f"Found {_n_zeroes} zero eigenvalues, expected {n_zeroes}"
                        "\n\nEnergies:\n" + vib_text
                    )
            else:
                raise ValueError(
                    f"Don't recognize number of imaginary modes: '{n_imaginary}'"
                )

            # Split the groups of vibrations
            real_energies = energies[n_zeroes + n_imaginary :]
            imaginary_energies = energies[:n_imaginary]
            zero_energies = energies[n_imaginary : n_zeroes + n_imaginary]
            is_imaginary = np.iscomplex(energies)

            imaginary_frequencies = np.abs(imaginary_energies) * Q_(1.0, "eV").m_as(
                "1/cm"
            )
            if n_imaginary > 0:
                text += f"There are {n_imaginary} imaginary frequencies:\n"
                frequencies = [f"  {v:9.1f}i" for v in imaginary_frequencies]
                text += "\n".join(frequencies)
                text += "\n\n"

            zero_frequencies = np.abs(zero_energies) * Q_(1.0, "eV").m_as("1/cm")
            text += f"There are {n_zeroes} frequencies that should be zero:\n"
            suffix = [
                "i" if v else ""
                for v in is_imaginary[n_imaginary : n_zeroes + n_imaginary]
            ]
            frequencies = [f"  {v:9.1f}{s}" for (v, s) in zip(zero_frequencies, suffix)]
            text += "\n".join(frequencies)
            text += "\nTheir magnitude gives a feel for the quality of the calculation."
            text += "\n\n"

            real_frequencies = np.abs(real_energies) * Q_(1.0, "eV").m_as("1/cm")
            n_max = _P["print frequencies"]
            if isinstance(n_max, str):
                if "no" in n_max.lower():
                    n_max = 0
                elif "all" in n_max.lower():
                    n_max = len(real_frequencies)
                else:
                    n_max = int(n_max)

            if n_max > 0 and len(real_frequencies) <= n_max:
                text += f"There are {len(real_frequencies)} true frequencies:\n"
                frequencies = [f"  {v:9.1f}" for v in real_frequencies]
                text += "\n".join(frequencies)
                text += "\n\n"

        if _P["spin multiplicity"] == "from system":
            spin = (self._working_configuration.spin_multiplicity - 1) / 2
        else:
            spin = (float(_P["spin multiplicity"]) - 1) / 2

        if _P["symmetry number"] == "from system":
            symmetry_number, axes = rotational_symmetry_number(
                atoms.get_chemical_symbols(), atoms.positions, atoms.get_masses()
            )
            if "linear" in axes:
                description = "linear, " + (
                    "centrosymmetric" if symmetry_number == 2 else "not centrosymmetric"
                )
            elif axes:
                description = ", ".join(
                    f"{count} C{order}" for order, count in axes.items()
                )
                description += " axes" if sum(axes.values()) > 1 else " axis"
            else:
                description = "no rotational symmetry"
            text += (
                f"The rotational symmetry number perceived from the structure is "
                f"{symmetry_number} ({description})."
            )
        else:
            symmetry_number = int(_P["symmetry number"])
            text += f"The rotational symmetry number is {symmetry_number}."
        text += "\n\n"

        # Store the symmetry number and frequencies
        data = {}
        if n_imaginary > 0:
            data["vibrational frequencies"] = [
                -v for v in imaginary_frequencies
            ] + list(real_frequencies)
        else:
            data["vibrational frequencies"] = list(real_frequencies)
        data["symmetry number"] = symmetry_number
        data["N saddle modes"] = n_imaginary
        if n_imaginary == 1:
            data["transition state frequency"] = imaginary_frequencies[0]

        # The thermodynamic functions, relative to the electronic energy (ASE
        # works in eV). The ZPE is included in ASE's U, H and G.
        thermo = IdealGasThermo(
            vib_energies=np.real(real_energies),
            potentialenergy=0.0,
            atoms=atoms,
            geometry=geometry,
            symmetrynumber=symmetry_number,
            spin=spin,
            vib_selection="all",
        )

        # The ZPE from the true vibrations only, excluding the nominally zero
        # translations and rotations and any imaginary modes
        ZPE = Q_(0.5 * float(np.sum(np.real(real_energies))), "eV").to("kJ/mol")
        printer.important(__(text, indent=4 * " "))
        text = f"Zero-point energy (ZPE) = {ZPE:.2f~#P} ({ZPE.to('kcal/mol'):.2f~#P})"
        ZPE = ZPE.magnitude
        data["ZPE"] = ZPE

        # The electronic energy of the undisplaced structure, from the
        # subflowchart's first calculation
        if self._E0 is not None:
            E0 = Q_(self._E0, "eV").m_as("kJ/mol")
            data["total_energy"] = Q_(E0, "kJ/mol").m_as("E_h")
            data["H 0"] = E0 + ZPE
            text += f"\nElectronic energy E = {E0:.3f} kJ/mol"
            text += f" ({data['total_energy']:.6f} E_h)"
            text += f"\nH(0 K) = E + ZPE = {data['H 0']:.3f} kJ/mol"
        else:
            E0 = None
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        Ts, Tunits = self._values_and_units(_P["T"], "K")
        Ps, Punits = self._values_and_units(_P["P"], "bar")

        # Conversion factors
        Efac = Q_(1, "eV").m_as("kJ/mol")
        Sfac = Q_(1, "eV/K").m_as("J/mol/K")
        Tfac = Q_(1, Tunits).m_as("K")
        Pfac = Q_(1, Punits).m_as("Pa")

        data["P"] = Ps
        data["T"] = Ts

        data["E thermal"] = []
        data["H thermal"] = []
        data["G thermal"] = []
        data["S"] = []
        if E0 is not None:
            data["U"] = []
            data["H"] = []
            data["G"] = []

        aP = []
        aT = []
        Ethermal = []
        Hthermal = []
        Gthermal = []
        Svalues = []
        for Porig in Ps:
            P = Pfac * Porig
            Etmp = []
            Htmp = []
            Gtmp = []
            Stmp = []
            for i, Torig in enumerate(Ts):
                T = Tfac * Torig
                aP.append(Porig if i == 0 else "")
                aT.append(Torig)
                U_value = Efac * thermo.get_internal_energy(T, verbose=False)
                H_value = Efac * thermo.get_enthalpy(T, verbose=False)
                S_value = Sfac * thermo.get_entropy(T, P, verbose=False)
                G_value = Efac * thermo.get_gibbs_energy(T, P, verbose=False)
                Etmp.append(U_value)
                Htmp.append(H_value)
                Gtmp.append(G_value)
                Stmp.append(S_value)
            data["E thermal"].append(Etmp)
            data["H thermal"].append(Htmp)
            data["G thermal"].append(Gtmp)
            data["S"].append(Stmp)
            if E0 is not None:
                data["U"].append([E0 + v for v in Etmp])
                data["H"].append([E0 + v for v in Htmp])
                data["G"].append([E0 + v for v in Gtmp])
            Ethermal.extend(Etmp)
            Hthermal.extend(Htmp)
            Gthermal.extend(Gtmp)
            Svalues.extend(Stmp)

        table = {}
        if len(Ps) == 1:
            title = (
                "Thermal corrections to the electronic energy (including the ZPE) at "
                f"{Ps[0]} {Punits}"
            )
        else:
            title = "Thermal corrections to the electronic energy (including the ZPE)"
            table[f"P ({Punits})"] = aP
        table[f"T ({Tunits})"] = aT
        table["U - E (kJ/mol)"] = Ethermal
        table["H - E (kJ/mol)"] = Hthermal
        table["S (J/mol/K)"] = Svalues
        table["G - E (kJ/mol)"] = Gthermal

        tmp = tabulate(
            table, headers="keys", tablefmt="rounded_outline", floatfmt=".3f"
        )

        length = len(tmp.splitlines()[0])
        text = "\n"
        text += title.center(length)
        text += "\n"
        text += tmp
        text += "\n"
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

        # Put any requested results into variables or tables
        self.store_results(
            data=data,
            create_tables=True,
            configuration=self._working_configuration,
        )

        # Citation!
        self.ase_read_bibliography()
        self.references.cite(
            raw=self._bibliography["ASE"],
            alias="ASE",
            module="structure_step",
            level=1,
            note="Main reference for ASE.",
        )

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "thermochemistry-step"
        parser = getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            # Any options for thermochemistry itself
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        # Now need to walk through the steps in the subflowchart...
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False, natoms=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        # Make sure the subflowchart has the data from the parent flowchart
        self.subflowchart.root_directory = self.flowchart.root_directory
        self.subflowchart.executor = self.flowchart.executor
        self.subflowchart.in_jobserver = self.subflowchart.in_jobserver

        if P is None:
            P = self.parameters.values_to_dict()

        result = self.header + "\n\n"

        text = (
            "Calculate the thermochemistry for the system in the harmonic approximation"
            " using the vibrational frequencies (phonons)."
        )
        if P.get("engine", "Model Chemistry") == "Model Chemistry":
            text += (
                " The energy and the Hessian will come from the flowchart's Model "
                "Chemistry: the analytic Hessian if the program provides one, "
                "otherwise finite differences of the gradients with a step size "
                "of {step size}."
            )
        else:
            text += (
                " The force constants will be calculated using the finite difference"
                " method, with a step size of {step size}, running the subflowchart"
                " for each displaced structure."
            )

        if P["# imaginary modes"] == "transition state":
            text += (
                " The system is a transition state, so the imaginary mode will be"
                " excluded from the calculation."
            )
        elif P["# imaginary modes"] == "perceive":
            text += " Any imaginary modes will be excluded from the calculation."
        else:
            text += (
                " {# imaginary modes} imaginary modes will be excluded from the "
                "calculation."
            )

        spin_multiplicity = P["spin multiplicity"]
        symmetry_number = P["symmetry number"]
        if spin_multiplicity == "from system":
            if symmetry_number == "from system":
                text += " The spin multiplicity will be taken from the system and the"
                text += " rotational symmetry number perceived from the structure."
            else:
                text += " The spin multiplicity will be taken from the system and the"
                text += f" symmetry number will be {symmetry_number}."
        else:
            if symmetry_number == "from system":
                text += " The symmetry number will be perceived from the structure and"
                text += f" the spin multiplicity will be {spin_multiplicity}."
            else:
                text += f" The spin multiplicity will be {spin_multiplicity} and the"
                text += f" symmetry number will be {symmetry_number}."

        if self.is_expr(P["print frequencies"]):
            text += " The frequencies will be printed depending on the value of"
            text += " {print frequencies}."
        elif P["print frequencies"] == "no":
            text += " The frequencies will not be printed."
        elif P["print frequencies"] == "all":
            text += " All the frequencies will be printed."
        else:
            text += " If there are no more than {print frequencies} frequencies,"
            text += " they will be printed."

        _T = P["T"]
        if isinstance(_T, tuple):
            _T = " ".join(_T)
        _P = P["P"]
        if isinstance(_P, tuple):
            _P = " ".join(_P)
        text += f" The thermodynamic properties will be calculated at {_T} and {_P}."

        result += str(__(text, indent=4 * " ", **P))

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        if not short and P.get("engine", "Model Chemistry") == "subflowchart":
            result += "\n\n"
            # Now walk through the steps in the subflowchart...
            while node is not None:
                try:
                    result += str(
                        __(node.description_text(), indent=7 * " ", wrap=False)
                    )
                except Exception as e:
                    print(f"Error describing structure flowchart: {e} in {node}")
                    self.logger.critical(
                        f"Error describing structure flowchart: {e} in {node}"
                    )
                    raise
                except:  # noqa: E722
                    print(
                        "Unexpected error describing structure flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    self.logger.critical(
                        "Unexpected error describing structure flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    raise
                result += "\n"
                node = node.next()

        return result

    def run(self):
        """Run a Thermochemistry step.

        Parameters
        ----------
        None

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)

        # Get the values of the parameters, dereferencing any variables
        _P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Have to fix formatting for printing...
        _PP = dict(_P)
        for key in _PP:
            if isinstance(_PP[key], units_class):
                _PP[key] = "{:~P}".format(_PP[key])

        # Get the  configuration
        _, self._working_configuration = self.get_system_configuration()
        n_atoms = self._working_configuration.n_atoms

        # Print what we are doing
        printer.important(
            __(
                self.description_text(_PP, short=True, natoms=n_atoms),
                indent=self.indent,
            )
        )

        step_size = _P["step size"].m_as("Å")

        # Where the energies and gradients come from. An old flowchart with a
        # subflowchart and no Model Chemistry keeps working.
        engine = _P.get("engine", "Model Chemistry")
        if engine == "Model Chemistry" and self.model_chemistry() is None:
            if self.subflowchart.get_node("1").next() is not None:
                printer.important(
                    __(
                        "There is no Model Chemistry before this step, so the "
                        "subflowchart is used for the energies and gradients.",
                        indent=self.indent + 4 * " ",
                    )
                )
                engine = "subflowchart"
            else:
                raise ValueError(
                    "Thermochemistry needs a Model Chemistry step before it, or a "
                    "subflowchart with the calculation of the energy and gradients."
                )

        tic = time.perf_counter_ns()
        if engine == "Model Chemistry":
            mc = self.model_chemistry()
            self.model = mc["level"]
            self._E0, self._vibrations, how = self.model_chemistry_vibrations(
                self._working_configuration, step_size
            )
            printer.important(
                __(
                    f"The Hessian was obtained with {mc['level']} ({how}).",
                    indent=self.indent + 4 * " ",
                )
            )
        else:
            # Get the Hessian ... save coordinates in case they are changed
            XYZ_save = self._working_configuration.atoms.coordinates
            self._XYZ0 = np.array(
                self._working_configuration.atoms.get_coordinates(fractionals=False)
            )
            self._E0 = None
            self._vibrations = self.run_ase_Hessian(
                step_size=step_size,
                on_error=_P["on error"],
                on_success=_P["on success"],
            )
            self._working_configuration.atoms.coordinates = XYZ_save
            # Name the properties after the engine's model chemistry
            self.model = self._subflowchart_model()
        toc = time.perf_counter_ns()
        self._results["t_elapsed"] = round((toc - tic) * 1.0e-9, 3)

        # Analyze the results
        self.analyze(_P)
        # Add other citations here or in the appropriate place in the code.
        # Add the bibtex to data/references.bib, and add a self.reference.cite
        # similar to the above to actually add the citation to the references.

        return next_node

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
