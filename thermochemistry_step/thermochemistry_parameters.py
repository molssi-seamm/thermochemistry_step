# -*- coding: utf-8 -*-
"""
Control parameters for the Thermochemistry step in a SEAMM flowchart
"""

import logging
import seamm

logger = logging.getLogger(__name__)


class ThermochemistryParameters(seamm.Parameters):
    """
    The control parameters for Thermochemistry.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"] : tuple
        A tuple of enumerated values.

    parameters["format_string"] : str
        A format string for "pretty" output.

    parameters["description"] : str
        A short string used as a prompt in the GUI.

    parameters["help_text"] : str
        A longer string to display as help for the user.

    See Also
    --------
    Thermochemistry, TkThermochemistry, Thermochemistry
    ThermochemistryParameters, ThermochemistryStep
    """

    parameters = {
        "engine": {
            "default": "Model Chemistry",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("Model Chemistry", "subflowchart"),
            "format_string": "",
            "description": "Energies and gradients from:",
            "help_text": (
                "The flowchart's Model Chemistry, run as a warm engine over MDI "
                "(the analytic Hessian when the program has one, else finite "
                "differences of the gradients), or this step's subflowchart, run "
                "once per displaced structure."
            ),
        },
        "approach": {
            "default": "Harmonic approximation",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "Harmonic approximation",
                "Isodesmic (CBH) with ATcT references",
            ),
            "format_string": "",
            "description": "Approach:",
            "help_text": "The approach or method for determining the thermochemistry.",
        },
        "cbh rungs": {
            "default": "all usable",
            "kind": "string",
            "default_units": "",
            "enumeration": ("all usable", "0", "1", "2", "3", "4", "1, 2", "2, 3"),
            "format_string": "",
            "description": "CBH rungs:",
            "help_text": (
                "The rungs of the connectivity-based hierarchy to use: every rung "
                "whose reference species are all in ATcT, or a list."
            ),
        },
        "protect": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "default",
                "none",
                "aromatic rings",
                "aromatic + small rings",
                "all rings",
            ),
            "format_string": "",
            "description": "Keep whole:",
            "help_text": (
                "Groups kept as single units when building the reactions. The "
                "default keeps aromatic ring systems, 3- and 4-membered rings and "
                "charge-separated or hypervalent groups (nitro, sulfonyl, ...)."
            ),
        },
        "cbh energies": {
            "default": "finite-difference harmonic",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "finite-difference harmonic",
                "existing results (run nothing)",
            ),
            "format_string": "",
            "description": "Species energies:",
            "help_text": (
                "Compute every species with the engine (optimization, then the "
                "Hessian for the harmonic thermochemistry), or read energies "
                "already stored as properties on the systems named CBH:<SMILES>."
            ),
        },
        "thermal corrections": {
            "default": "harmonic",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "harmonic",
                "none (the energy is the enthalpy, e.g. MOPAC)",
            ),
            "format_string": "",
            "description": "Thermal corrections:",
            "help_text": (
                "Harmonic thermal corrections, or none when the engine's energy is "
                "already an enthalpy of formation at 298.15 K (MOPAC and other "
                "NDDO methods)."
            ),
        },
        "initial conformer": {
            "default": "lowest of the embeddings (MMFF)",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("lowest of the embeddings (MMFF)", "single embedding"),
            "format_string": "",
            "description": "Fragments start from:",
            "help_text": (
                "The starting structure of each fragment built from its SMILES: the "
                "lowest MMFF94 conformer of several RDKit embeddings (all-anti chains, "
                "as CBH-ANL starts from the lowest rotamer), or a single embedding."
            ),
        },
        "conformer embeddings": {
            "default": 20,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of embeddings:",
            "help_text": "How many RDKit embeddings to rank by MMFF94.",
        },
        "target conformer": {
            "default": "as given",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("as given", "lowest of the embeddings (MMFF)"),
            "format_string": "",
            "description": "Target starts from:",
            "help_text": (
                "Keep the target structure as given, or rebuild it from its SMILES "
                "as the lowest MMFF94 conformer of several embeddings before the "
                "optimization (a single embedding of a long chain is often gauche)."
            ),
        },
        "optimize fragments": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Optimize the structures:",
            "help_text": (
                "Whether to optimize each species (ASE BFGS over the engine) "
                "before the Hessian."
            ),
        },
        "fragment convergence": {
            "default": 0.005,
            "kind": "float",
            "default_units": "eV/Å",
            "enumeration": tuple(),
            "format_string": ".3f",
            "description": "Optimization convergence:",
            "help_text": "The maximum force for a converged optimization.",
        },
        "fragment max steps": {
            "default": 200,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Maximum optimization steps:",
            "help_text": "The maximum number of optimization steps per species.",
        },
        "energy property": {
            "default": "total energy#*#*",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Energy property:",
            "help_text": (
                "The stored property holding each species' energy (a glob pattern "
                "matching exactly one property)."
            ),
        },
        "ZPE property": {
            "default": "ZPE#*#*",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "ZPE property:",
            "help_text": "The stored property holding the zero-point energy.",
        },
        "enthalpy property": {
            "default": "H(298.15 K)#*#*",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Enthalpy property:",
            "help_text": (
                "The stored property holding the absolute enthalpy at 298.15 K "
                "(its level of theory may differ from the energy's)."
            ),
        },
        "low-level energy property": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Energy at the enthalpy's level:",
            "help_text": (
                "If the enthalpy comes from a different level of theory than the "
                "energy, the energy at that level, so that H - E can be formed. "
                "Blank: the same as the energy property."
            ),
        },
        "fragment store": {
            "default": "use if available",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("use if available", "ignore"),
            "format_string": "",
            "description": "Fragment store:",
            "help_text": (
                "Use the energies of species already in the fragment-energy store "
                "(seamm-thermochemistry-fragments) at this model chemistry instead "
                "of computing them, or ignore the store."
            ),
        },
        "fragment store file": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "",
            "description": "Fragment store file:",
            "help_text": ("The store: fragments.db next to the ATcT cache by default."),
        },
        "reference values": {
            "default": "ATcT",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "ATcT",
                "laddered",
                "ATcT below the threshold, else laddered",
            ),
            "format_string": "",
            "description": "Reference values:",
            "help_text": (
                "Where the reference species' enthalpies of formation come from: "
                "ATcT for every species; a laddered reference set (CBH-ANL: ATcT "
                "only for the CBH-0 species, the rest from their own lower-rung "
                "reactions at a high level, built with seamm-thermochemistry-ladder); "
                "or ATcT where its 95 % uncertainty is below the threshold and the "
                "ladder otherwise."
            ),
        },
        "ladder": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Ladder:",
            "help_text": "The name of the laddered reference set.",
        },
        "ATcT uncertainty threshold": {
            "default": 0.5,
            "kind": "float",
            "default_units": "kJ/mol",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "ATcT uncertainty threshold:",
            "help_text": (
                "Use the ATcT value when its 95 % uncertainty is at or below this."
            ),
        },
        "ATcT version": {
            "default": "latest",
            "kind": "string",
            "default_units": "",
            "enumeration": ("latest",),
            "format_string": "",
            "description": "ATcT version:",
            "help_text": "The version of the Thermochemical Network to use.",
        },
        "ATcT directory": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "",
            "description": "ATcT directory:",
            "help_text": (
                "Where the ATcT versions are kept: the installation's "
                "Parameters/thermochemistry/ATcT by default."
            ),
        },
        "covariances": {
            "default": "cached only",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("cached only", "fetch if needed"),
            "format_string": "",
            "description": "ATcT covariances:",
            "help_text": (
                "Use only the covariance blocks already cached, or fetch missing "
                "ones from atct.anl.gov (not on a compute node without network)."
            ),
        },
        "compare with ATcT": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Compare with ATcT:",
            "help_text": "Report the difference if the target itself is in ATcT.",
        },
        "hessian route": {
            "default": "automatic",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "automatic",
                "warm engine over MDI",
                "separate tasks",
            ),
            "format_string": "",
            "description": "Energies and Hessian via:",
            "help_text": (
                "How the Model Chemistry is driven: automatically (the engine when the "
                "program has an analytic Hessian or no task route, else finite "
                "differences as separate tasks where the job runs them); the warm MDI "
                "engine, which also allows the optimization of the structure; or "
                "separate tasks. The CBH approach uses the engine by default, since it "
                "optimizes each species."
            ),
        },
        "step size": {
            "default": 0.01,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Step size:",
            "help_text": (
                "The size of the step for finite differences (not used when the "
                "program provides the analytic Hessian)."
            ),
        },
        "# imaginary modes": {
            "default": "perceive",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("perceive", "transition state"),
            "format_string": "",
            "description": "Number of imaginary modes:",
            "help_text": "The number of imaginary modes to ignore.",
        },
        "print frequencies": {
            "default": "all",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("all", "none"),
            "format_string": "",
            "description": "Maximum frequencies to print:",
            "help_text": "Whether to print the frequencies of vibration.",
        },
        "T": {
            "default": "200, 250, 298.15, 300:1000:50",
            "kind": "string",
            "default_units": "K",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Temperatures:",
            "help_text": "The temperatures for the thermodynamic functions.",
        },
        "P": {
            "default": "1.0",
            "kind": "string",
            "default_units": "bar",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Pressures:",
            "help_text": "The pressures for the thermodynamic functions.",
        },
        "spin multiplicity": {
            "default": "from system",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "from system",
                "singlet",
                "doublet",
                "triplet",
                "quartet",
                "quintet",
                "sextet",
                "septet",
                "octet",
            ),
            "format_string": "",
            "description": "Spin multiplicity:",
            "help_text": "The spin multiplicity of the system as word or integer.",
        },
        "symmetry number": {
            "default": "from system",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("from system",),
            "format_string": "",
            "description": "Symmetry number:",
            "help_text": "The symmetry number of the system.",
        },
        "on success": {
            "default": "keep last subdirectory",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "keep last subdirectory",
                "keep all subdirectories",
                "delete all subdirectories",
            ),
            "format_string": "",
            "description": "On success:",
            "help_text": "Which subdirectories to keep.",
        },
        "on error": {
            "default": "keep all subdirectories",
            "kind": "enum",
            "default_units": "",
            "enumeration": (
                "keep last subdirectory",
                "keep all subdirectories",
                "delete all subdirectories",
            ),
            "format_string": "",
            "description": "On error:",
            "help_text": "Which subdirectories to keep if there is an error.",
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": None,
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": "The results to save to variables or in tables.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ThermochemistryParameters.__init__")

        super().__init__(
            defaults={
                **ThermochemistryParameters.parameters,
                **defaults,
            },
            data=data,
        )
