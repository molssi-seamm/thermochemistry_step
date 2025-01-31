# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from Thermochemistry
"""

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
        "dimensionality": "scalar",
        "property": "U#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "H": {
        "description": "The enthalpy of the system",
        "dimensionality": "scalar",
        "property": "H#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "S": {
        "description": "The entropy of the system",
        "dimensionality": "scalar",
        "property": "S#Thermochemistry#{model}",
        "type": "float",
        "units": "kJ/mol/K",
        "format": ".6f",
    },
    "G": {
        "description": "The free energy of the system",
        "dimensionality": "scalar",
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
    "symmetry number": {
        "description": "Symmetry number",
        "dimensionality": "scalar",
        "type": "integer",
        "format": "",
    },
    "E thermal": {
        "description": "the thermal correction to the electronic energy",
        "dimensionality": "[nPs*[nTs]]",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "H thermal": {
        "description": "the enthalpy correction to the electronic energy",
        "dimensionality": "[nPs*[nTs]]",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
    "G thermal": {
        "description": "the free energy correction to the electronic energy",
        "dimensionality": "[nPs*[nTs]]",
        "type": "float",
        "units": "kJ/mol",
        "format": ".6f",
    },
}
