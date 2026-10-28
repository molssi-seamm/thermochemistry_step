.. _user-guide:

**********
User Guide
**********

The Thermochemistry step calculates the vibrational frequencies of a molecule by
finite differences of the gradients, and from them the ideal-gas, rigid-rotor,
harmonic-oscillator thermodynamic functions at the requested temperatures and
pressures.

How it works
============

The step holds a subflowchart with the engine that provides the energy and
gradients: a MOPAC, ORCA, Gaussian, Psi4, LAMMPS or MLFF energy step, for example.
The subflowchart runs once for the undisplaced structure and once for each
displacement of each atom in ±x, ±y and ±z (6N+1 runs). The engine step must pass
its energy and gradients back as JSON: in its *Results* tab, store ``energy`` and
``gradients`` as JSON. Without these the step stops with "There are no energy and
gradients in properties.json".

The frequencies come from ASE's finite-difference Hessian. The six (five for a
linear molecule) nominally zero translations and rotations are reported so that
the quality of the gradients can be judged, and imaginary modes are excluded from
the thermodynamic functions: by default the step perceives how many there are, or
you can declare the structure a transition state or give the number.

Parameters
==========

Step size
    The finite-difference displacement, 0.01 Å by default.

Number of imaginary modes
    ``perceive`` (default), ``transition state`` (one), or an integer.

Temperatures and pressures
    Lists in the usual SEAMM form, e.g. ``200, 250, 298.15, 300:1000:50`` K and
    ``1.0`` bar.

Spin multiplicity
    Taken from the system by default. It contributes R ln(2S+1) to the entropy.

Symmetry number
    The rotational symmetry number σ, which enters the entropy as -R ln σ. By
    default (``from system``) it is perceived from the structure by testing
    rotation axes against the geometry; the step prints the number and the axes
    it found, e.g. "12 (4 C3, 3 C2 axes)" for methane. A geometry must be
    symmetric to within 0.1 Å for an axis to be found, which optimized structures
    are. Give the number explicitly if the structure is deliberately distorted.

Results
=======

The step prints the ZPE, the electronic energy E of the undisplaced structure,
H(0 K) = E + ZPE, and a table of the thermal corrections U-E, H-E and G-E (each
including the ZPE, in kJ/mol) and the entropy S (J/mol/K) at every temperature
and pressure.

The results that can be stored in variables, tables or as properties of the
configuration are the frequencies, the symmetry number, the number of saddle
modes and the transition-state frequency, the ZPE, the total (electronic) energy,
H(0 K), and U, H, S and G. The last four are lists over the pressures and then
the temperatures, in the same order as the stored ``P`` and ``T`` lists. The
properties are named with the engine's model chemistry, for example
``H#Thermochemistry#PM6-ORG`` or ``S#Thermochemistry#B3LYP/def2-TZVP``.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
