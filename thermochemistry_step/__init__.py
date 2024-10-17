# -*- coding: utf-8 -*-

"""
thermochemistry_step
A SEAMM plug-in for Thermochemistry
"""

# Bring up the classes so that they appear to be directly in
# the thermochemistry_step package.

from .thermochemistry import Thermochemistry  # noqa: F401, E501
from .thermochemistry_step import ThermochemistryStep  # noqa: F401, E501
from .tk_thermochemistry import TkThermochemistry  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
