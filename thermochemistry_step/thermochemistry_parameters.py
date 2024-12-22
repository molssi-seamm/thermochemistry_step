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
        "approach": {
            "default": "Harmonic approximation",
            "kind": "enum",
            "default_units": "",
            "enumeration": ("Harmonic approximation",),
            "format_string": "",
            "description": "Approach:",
            "help_text": "The approach or method for determining the thermochemistry.",
        },
        "step size": {
            "default": 0.01,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".g",
            "description": "Step size:",
            "help_text": "The size of the step for finite differences.",
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
