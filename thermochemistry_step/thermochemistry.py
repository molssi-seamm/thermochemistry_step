# -*- coding: utf-8 -*-

"""Non-graphical part of the Thermochemistry step in a SEAMM flowchart
"""

import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401
import sys
import time

from ase.thermochemistry import IdealGasThermo
import numpy as np
from tabulate import tabulate

import thermochemistry_step
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
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class Thermochemistry(seamm.Node, ASE_mixin):
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

        if n_atoms == 1:
            symmetry_number = 1
            real_frequencies = []
        else:
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
                    n_max = len(frequencies)

            if n_max > 0 and len(frequencies) <= n_max:
                text += f"There are {len(real_frequencies)} true frequencies:\n"
                frequencies = [f"  {v:9.1f}" for v in real_frequencies]
                text += "\n".join(frequencies)
                text += "\n\n"

        potential_energy = 0.0
        if _P["spin multiplicity"] == "from system":
            spin = (self._working_configuration.spin_multiplicity - 1) / 2
        else:
            spin = (float(_P["spin multiplicity"]) - 1) / 2

        if _P["symmetry number"] == "from system":
            raise NotImplementedError("symmetry number from system")
        else:
            symmetry_number = int(_P["symmetry number"])

        thermo = IdealGasThermo(
            vib_energies=energies.real,
            potentialenergy=potential_energy,
            atoms=atoms,
            geometry=geometry,
            symmetrynumber=symmetry_number,
            spin=spin,
        )
        # If have imaginary frequencies, patch up those in ASE
        thermo.vib_energies = np.real(real_energies)

        ZPE = Q_(self._vibrations.get_zero_point_energy(), "eV").to("kJ/mol")
        text += f"Zero-point energy (ZPE) = {ZPE:.1f~#P} ({ZPE.to('kcal/mol'):.1f~#P})"
        printer.important(__(text, indent=4 * " "))

        H = []
        S = []
        G = []

        T_list = _P["T"]
        P_list = _P["P"]

        if isinstance(T_list, tuple):
            # (list, units)
            Ts, Tunits = T_list
            Ts = parse_list(Ts)
        else:
            Tunits = str(Ts.u)
            Ts = [T_list.magnitude]
        if isinstance(P_list, tuple):
            # (list, units)
            Ps, Punits = P_list
            Ps = parse_list(Ps)
        else:
            Ps = [P_list.magnitude]
            Punits = str(P_list.u)

        # Conversion factor
        Efac = Q_(1, "eV").m_as("kcal/mol")
        Sfac = Q_(1, "eV/K").m_as("cal/mol/K")
        Tfac = Q_(1, Tunits).m_as("K")
        Pfac = Q_(1, Punits).m_as("Pa")

        if len(Ps) == 1:
            title = (
                f"ZPE and thermal contributions to thermodynamic functions at {Ps[0]}"
                f" {Punits}"
            )
            P = Pfac * Ps[0]
            for T in Ts:
                H.append(Efac * thermo.get_enthalpy(T, verbose=False))
                S.append(Sfac * thermo.get_entropy(T, P, verbose=False))
                G.append(Efac * thermo.get_gibbs_energy(T, P, verbose=False))

            table = {}
            table[f"T ({Tunits})"] = Ts
            table["H (kcal/mol)"] = H
            table["S (cal/mol/K)"] = S
            table["G (kcal/mol)"] = G
        elif len(Ts) == 1:
            title = (
                f"ZPE and thermal contributions to thermodynamic functions at {Ts[0]}"
                f" {Tunits}"
            )
            T = Tfac * Ts[0]
            for P in Ps:
                H.append(Efac * thermo.get_enthalpy(T, verbose=False))
                S.append(Sfac * thermo.get_entropy(T, P, verbose=False))
                G.append(Efac * thermo.get_gibbs_energy(T, P, verbose=False))

            table = {}
            table[f"P ({Punits})"] = Ps
            table["H (kcal/mol)"] = H
            table["S (cal/mol/K)"] = S
            table["G (kcal/mol)"] = G
        else:
            title = "ZPE and thermal contributions to thermodynamic functions"
            aP = []
            aT = []
            for Porig in Ps:
                P = Pfac * Porig
                for i, T in enumerate(Ts):
                    if i == 0:
                        aP.append(Porig)
                    else:
                        aP.append("")
                    aT.append(T)
                    T *= Tfac
                    H.append(Efac * thermo.get_enthalpy(T, verbose=False))
                    S.append(Sfac * thermo.get_entropy(T, P, verbose=False))
                    G.append(Efac * thermo.get_gibbs_energy(T, P, verbose=False))

            table = {}
            table[f"P ({Punits})"] = aP
            table[f"T ({Tunits})"] = aT
            table["H (kcal/mol)"] = H
            table["S (cal/mol/K)"] = S
            table["G (kcal/mol)"] = G

        tmp = tabulate(
            table, headers="keys", tablefmt="rounded_outline", floatfmt=".2f"
        )

        length = len(tmp.splitlines()[0])
        text = "\n"
        text += title.center(length)
        text += "\n"
        text += tmp
        text += "\n"
        printer.important(__(text, indent=4 * " ", wrap=False, dedent=False))

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
            " using the vibrational frequencies (phonons). The force constants will be"
            " calculated using the finite difference method, with a step size of"
            " {step size}."
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
                text += " The spin multiplicity and symmetry number will be taken from"
                text += " the system."
            else:
                text += " The spin multiplicity will be taken from the system and the"
                text += f" symmetry number will be {symmetry_number}."
        else:
            if symmetry_number == "from system":
                text += " The symmetry number will be taken from the system and the"
                text += f" spin multiplicity will be {spin_multiplicity}."
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

        if not short:
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

        # Get the Hessian ... save coordinates in case they are changed
        XYZ_save = self._working_configuration.atoms.coordinates
        tic = time.perf_counter_ns()
        self._vibrations = self.run_ase_Hessian(step_size=step_size)
        toc = time.perf_counter_ns()
        self._results["t_elapsed"] = round((toc - tic) * 1.0e-9, 3)
        self._working_configuration.atoms.coordinates = XYZ_save

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
