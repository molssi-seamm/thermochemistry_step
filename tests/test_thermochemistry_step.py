#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `thermochemistry_step` package."""

import pytest  # noqa: F401
import thermochemistry_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = thermochemistry_step.Thermochemistry()
    assert (
        str(type(result))
        == "<class 'thermochemistry_step.thermochemistry.Thermochemistry'>"
    )
