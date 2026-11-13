# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from Thermochemistry"""

metadata = {}

"""Description of the computational models for Thermochemistry.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the Thermochemistry keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that Thermochemistry produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    "total_energy": {
        "description": "The total energy",
        "dimensionality": "scalar",
        "property": "total energy#Thermochemistry#{model}",
        "type": "float",
        "units": "E_h",
    },
    "H 0": {
        "description": "the electronic energy plus ZPE",
        "dimensionality": "scalar",
        "property": "H(0 K)#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "U": {
        "description": "the internal energy",
        "dimensionality": "[nPs*[nTs]]",
        "property": "U#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "H": {
        "description": "The enthalpy of the system",
        "dimensionality": "[nPs*[nTs]]",
        "property": "H#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "S": {
        "description": "The entropy of the system",
        "dimensionality": "[nPs*[nTs]]",
        "property": "S#Thermochemistry#{model}",
        "type": "float",
        "units": "J/mol/K",
        "format": ".6f",
    },
    "G": {
        "description": "The free energy of the system",
        "dimensionality": "[nPs*[nTs]]",
        "property": "G#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "ZPE": {
        "description": "the zero point energy",
        "dimensionality": "scalar",
        "property": "ZPE#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "ZPE uncorrected": {
        "description": "the zero point energy from the uncorrected frequencies",
        "dimensionality": "scalar",
        "property": "ZPE uncorrected#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "frequency correction status": {
        "description": "whether the frequency correction applied to every mode",
        "dimensionality": "scalar",
        "property": "frequency correction status#Thermochemistry#{model}",
        "type": "string",
        "units": "",
        "format": "s",
    },
    "frequency correction": {
        "description": "the correction of the frequencies (id, version and range)",
        "dimensionality": "scalar",
        "property": "frequency correction#Thermochemistry#{model}",
        "type": "string",
        "units": "",
        "format": "s",
    },
    "P": {
        "description": "Pressure",
        "dimensionality": "[nPs]",
        "property": "P#Thermochemistry#{model}",
        "type": "float",
        "units": "atm",
        "format": ".2f",
    },
    "T": {
        "description": "Temperature",
        "dimensionality": "[nTs]",
        "property": "T#Thermochemistry#{model}",
        "type": "float",
        "units": "K",
        "format": ".2f",
    },
    "DfH 0K": {
        "description": "The enthalpy of formation at 0 K (CBH with ATcT references)",
        "dimensionality": "scalar",
        "property": "DfH(0 K)#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "DfH 298K": {
        "description": (
            "The enthalpy of formation at 298.15 K (CBH with ATcT references)"
        ),
        "dimensionality": "scalar",
        "property": "DfH(298.15 K)#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "DfH uncertainty": {
        "description": (
            "The 95 % uncertainty of the ATcT reference part of the enthalpy of "
            "formation"
        ),
        "dimensionality": "scalar",
        "property": "DfH uncertainty (95%)#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "DfH uncertainty method": {
        "description": (
            "How the uncertainty was obtained: from the ATcT covariances, or in "
            "quadrature without them, plus any ladder method term"
        ),
        "dimensionality": "scalar",
        "property": "DfH uncertainty method#Thermochemistry#{model}",
        "type": "string",
        "format": "",
    },
    "library DfH 298K": {
        "description": "The target's own value at 298.15 K in the reference library",
        "dimensionality": "scalar",
        "property": "library DfH 298K#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "library DfH 0K": {
        "description": "The target's own value at 0 K in the reference library",
        "dimensionality": "scalar",
        "property": "library DfH 0K#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "library DfH uncertainty": {
        "description": (
            "The 95 % uncertainty of the target's library value, with its "
            "source's offset"
        ),
        "dimensionality": "scalar",
        "property": "library DfH uncertainty#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "library source": {
        "description": "The source of the target's library value",
        "dimensionality": "scalar",
        "property": "library source#Thermochemistry#{model}",
        "type": "string",
        "format": "",
    },
    "CBH - library": {
        "description": (
            "The CBH enthalpy of formation (highest rung, 298.15 K, else "
            "0 K) minus the target's library value"
        ),
        "dimensionality": "scalar",
        "property": "CBH - library#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "library check": {
        "description": (
            "Whether the CBH value agrees with the target's library value "
            "within twice the library's uncertainty"
        ),
        "dimensionality": "scalar",
        "property": "library check#Thermochemistry#{model}",
        "type": "string",
        "format": "",
    },
    "CBH references": {
        "description": (
            "Where the reference values of the reaction used came from: ATcT-only, "
            "laddered, mixed (the reference library too), mixed: several sources "
            "(for one element), computed references or user references"
        ),
        "dimensionality": "scalar",
        "property": "CBH references#Thermochemistry#{model}",
        "type": "string",
        "format": "",
    },
    "references": {
        "description": (
            "The versions of ATcT and the reference libraries, and the library "
            "values used (id, source, citation, conversion), for reproducing "
            "and auditing the result"
        ),
        "dimensionality": "dictionary",
        "type": "json",
        "format": "",
    },
    "CBH rung": {
        "description": "The CBH rung the enthalpy of formation comes from",
        "dimensionality": "scalar",
        "property": "CBH rung#Thermochemistry#{model}",
        "type": "integer",
        "format": "",
    },
    "CBH reaction": {
        "description": "The CBH reaction used",
        "dimensionality": "scalar",
        "property": "CBH reaction#Thermochemistry#{model}",
        "type": "string",
        "format": "",
    },
    "CBH energy level": {
        "description": "The model chemistry of the species' energies",
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "CBH geometry level": {
        "description": (
            "The model chemistry of the geometries, frequencies, ZPE and thermal "
            "corrections"
        ),
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "CBH high level": {
        "description": (
            "How the 'high' preset chose the energy level for this target "
            "(canonical or DLPNO CCSD(T)-F12, or why there is no 'high' result)"
        ),
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "CBH high canonical limit": {
        "description": (
            "The largest canonical CCSD(T)-F12 the 'high' preset allowed (cc-pVTZ-F12 "
            "basis functions)"
        ),
        "dimensionality": "scalar",
        "type": "integer",
        "format": "",
    },
    "CBH composite recipe": {
        "description": (
            "The composite method (SEAMM:CMP@<recipe>) whose geometry and energy "
            "levels were used, if any"
        ),
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "CBH composite recipe version": {
        "description": "The version of the composite recipe",
        "dimensionality": "scalar",
        "type": "integer",
        "format": "",
    },
    "CBH composite recipe fingerprint": {
        "description": "The fingerprint of the composite recipe (sha256)",
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "ZPE scale factor": {
        "description": "The factor the harmonic ZPEs were scaled by",
        "dimensionality": "scalar",
        "type": "float",
        "format": ".4f",
    },
    "ZPE scale factor source": {
        "description": "Where the ZPE scale factor comes from, and any caution",
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "CBH rung spread": {
        "description": (
            "The range of the enthalpy of formation at 298.15 K over the CBH rungs "
            "from CBH-1 up, an indicator of its uncertainty"
        ),
        "dimensionality": "scalar",
        "property": "CBH rung spread#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "ATcT version": {
        "description": "The ATcT Thermochemical Network version used",
        "dimensionality": "scalar",
        "type": "string",
        "format": "",
    },
    "ATcT DfH 0K": {
        "description": "ATcT's own enthalpy of formation of the target at 0 K",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "ATcT DfH 298K": {
        "description": "ATcT's own enthalpy of formation of the target at 298.15 K",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "ATcT DfH uncertainty": {
        "description": "ATcT's 95 % uncertainty for the target at 298.15 K",
        "dimensionality": "scalar",
        "type": "float",
        "units": "kJ/mol",
        "format": ".2f",
    },
    "CBH results": {
        "description": "The results for every usable rung",
        "dimensionality": "dictionary",
        "type": "json",
        "format": "",
    },
    "transition state frequency": {
        "description": "Transition state curvature",
        "dimensionality": "scalar",
        "property": "transition state frequency#Thermochemistry#{model}",
        "type": "float",
        "units": "cm^-1",
        "format": ".1f",
    },
    "N saddle modes": {
        "description": "Number of saddle point modes",
        "dimensionality": "scalar",
        "property": "number of saddle modes#Thermochemistry#{model}",
        "type": "integer",
        "format": "",
    },
    "vibrational frequencies": {
        "description": "Vibrational frequencies",
        "dimensionality": "[nvibs]",
        "property": "vibrational frequencies#Thermochemistry#{model}",
        "type": "float",
        "units": "cm^-1",
        "format": ".1f",
    },
    "corrected vibrational frequencies": {
        "description": "Vibrational frequencies after the frequency correction",
        "dimensionality": "[nvibs]",
        "property": "corrected vibrational frequencies#Thermochemistry#{model}",
        "type": "float",
        "units": "cm^-1",
        "format": ".1f",
    },
    "symmetry number": {
        "description": "Symmetry number",
        "dimensionality": "scalar",
        "type": "integer",
        "format": "",
    },
    "E thermal": {
        "description": (
            "U(T) - E, the thermal correction to the electronic energy, including "
            "the ZPE"
        ),
        "dimensionality": "[nPs*[nTs]]",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "H thermal": {
        "description": (
            "H(T) - E, the enthalpy correction to the electronic energy, including "
            "the ZPE"
        ),
        "dimensionality": "[nPs*[nTs]]",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "G thermal": {
        "description": (
            "G(T) - E, the free energy correction to the electronic energy, "
            "including the ZPE"
        ),
        "dimensionality": "[nPs*[nTs]]",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
}
