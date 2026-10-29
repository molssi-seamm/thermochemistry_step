.. _user-guide:

**********
User Guide
**********

The Thermochemistry step has two approaches. The **harmonic approximation**
calculates the vibrational frequencies of a molecule from its Hessian, and from
them the ideal-gas, rigid-rotor, harmonic-oscillator thermodynamic functions at
the requested temperatures and pressures. The **isodesmic (CBH) approach with
ATcT references** gives the enthalpy of formation of a molecule at 0 K and
298.15 K from the connectivity-based hierarchy of error-cancelling reactions,
with the reference species taken from the Active Thermochemical Tables.

How it works
============

The step needs the energy and the Hessian (second derivatives) of the
structure. By default they come from the flowchart's **Model Chemistry**: put a
Model Chemistry step before this one to choose the program and method (MOPAC
PM6-ORG, xTB, an ORCA DFT model, an MLFF model, ...). The program runs as a
warm engine over MDI. When it provides an analytic Hessian, as ORCA does for
Hartree-Fock and most density functionals, that is used; otherwise the Hessian
is built from central finite differences of the gradients over the same
connection, or as separate calculations where the job runs its tasks. The step
reports which route it took. The engine is started with the environment the
program's own runs get: for a program installed with environment modules
(``installation = modules`` in its ``.ini`` file, as on many clusters) its
modules are loaded first.

The alternative, "Energies and gradients from: subflowchart", keeps a
subflowchart inside the step holding the engine (a MOPAC, ORCA, Gaussian, Psi4,
LAMMPS or MLFF energy step, for example). The subflowchart runs once for the
undisplaced structure and once for each displacement of each atom in ±x, ±y and
±z (6N+1 runs). The engine step must pass its energy and gradients back as JSON:
in its *Results* tab, store ``energy`` and ``gradients`` as JSON. Without these
the step stops with "There are no energy and gradients in properties.json". Use
this route for calculations an MDI engine cannot do. A flowchart made before
the Model Chemistry route existed, with a subflowchart and no Model Chemistry,
still runs as before.

The frequencies come from the Hessian through ASE. The six (five for a linear
molecule) nominally zero translations and rotations are reported so that the
quality of the gradients can be judged, and imaginary modes are excluded from
the thermodynamic functions: by default the step perceives how many there are,
or you can declare the structure a transition state or give the number.

Parameters
==========

Energies and gradients from
    The flowchart's Model Chemistry (default) or this step's subflowchart.

Energies and Hessian via
    ``automatic`` (default for the harmonic approach): the warm MDI engine when
    the program has an analytic Hessian or no task route, else finite
    differences as separate tasks where the job runs them; ``warm engine over
    MDI``; or ``separate tasks``. The CBH approach uses the engine by default,
    since it optimizes each species over it.

Step size
    The finite-difference displacement, 0.01 Å by default; not used when the
    program provides the analytic Hessian.

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

Enthalpies of formation: the CBH approach
=========================================

Choose "Approach: Isodesmic (CBH) with ATcT references". The step builds the
connectivity-based hierarchy (CBH) reactions of the molecule from its bonds, so
the structure must have bonds (built from SMILES, or read from an SDF or MOL
file). CBH-0 breaks the molecule into hydrides of its atoms, CBH-1 into its
bonds, CBH-2 keeps each atom with all its neighbours, CBH-3 each bond with its
neighbours and CBH-4 each atom with two shells of neighbours: the higher the
rung, the more of the chemical environment, and of the error of the method, is
the same on both sides of the reaction. The reaction enthalpy is calculated with
the engine for the target and every reference species, and the enthalpy of
formation of the target follows from it and the ATcT values of the reference
species, with the 95 % uncertainty of the ATcT part (using the ATcT covariances
when they are cached). If the target itself is in ATcT the difference is
reported too.

Radicals and ions are handled: every fragment containing the radical or
charged atom keeps that spin or charge, and the ATcT entry is selected with
the same charge and multiplicity. The reactions follow from the SMILES of the
structure, so the system's charge and spin multiplicity must agree with it;
the step stops if they do not. Build a radical from its SMILES
(``[CH2]C`` for the ethyl radical), which sets both.

The ATcT tables are not distributed with SEAMM: install a version first with
``seamm-thermochemistry-atct download`` (see the seamm_thermochemistry
documentation). A rung is usable only when all its reference species are in
ATcT, and the step reports the rungs it could not use and why; the result
stored is the highest usable rung, and every rung goes into the table ``CBH``.

CBH rungs
    ``all usable`` (default), or a list such as ``1, 2``.

Keep whole
    Groups kept as single units when the reactions are built, so that they are
    never broken: by default aromatic ring systems, 3- and 4-membered rings and
    charge-separated or hypervalent groups (nitro, sulfonyl, phosphoryl, PF6-,
    ...). ``none`` builds the plain hierarchy.

Species energies
    ``finite-difference harmonic`` (default): each species is built from its
    SMILES as a system named ``CBH:<SMILES>``, optimized with the engine and
    its Hessian calculated, giving H(0 K) = E + ZPE and H(298.15 K). The
    optimization convergence and maximum steps, and the finite-difference step
    size, are as for the harmonic approach. Or ``existing results (run
    nothing)``: the energies are read from properties already stored on those
    systems, named by the three property patterns (energy, ZPE and enthalpy;
    each glob must match exactly one property), for example from an earlier
    job that computed the fragments at a higher level.

Fragments start from, Target starts from
    Each fragment built from its SMILES starts, by default, from the lowest
    MMFF94 conformer of 20 RDKit embeddings (an n-alkane all-anti, ethanol
    anti), as CBH-ANL starts every species from its lowest rotamer; a single
    embedding of a long chain is often gauche. The target is kept as given
    unless told to be rebuilt the same way. The chosen chain dihedrals are
    printed.

Thermal corrections
    ``harmonic`` (default), or ``none`` when the engine's energy is already an
    enthalpy of formation at 298.15 K, as for MOPAC and other NDDO methods:
    then no Hessian is calculated and only the 298.15 K result is available.
    Note that with NDDO methods CBH-0 is of no use, since their H2 is far from
    experiment.

Reference values
    ``ATcT`` for every reference species; ``laddered``: a CBH-ANL reference
    set built with ``seamm-thermochemistry-ladder`` (ATcT only for the CBH-0
    species, every other species from its own lower-rung reaction at a high
    level), named in *Ladder*; or ``ATcT below the threshold, else laddered``.
    With a ladder the uncertainty is the ATcT part through the ladder's
    coefficients on the anchors (cᵀΣc with the covariance) plus the ladder's
    measured method term, reported separately. Validate laddered results at
    0 K: ATcT's 298.15 K values are conformer averages.

Fragment store
    Species already in the fragment-energy store at this model chemistry
    (``seamm-thermochemistry-fragments harvest`` promotes them there from
    finished jobs) are taken from it instead of being computed, and serve as
    the fallback in "existing results" mode. The species table says where
    each value came from; a record without a settings fingerprint is
    "unverified" and the output says so. When the engine provides them (ORCA),
    the settings fingerprints of each species' calculations are stored with its
    energies, so the harvest files verified records.

ATcT version, directory, covariances
    The version of the Thermochemical Network (``latest`` cached by default),
    where the versions are kept, and whether to fetch missing covariance blocks
    from atct.anl.gov (``cached only`` by default: a compute node usually has
    no network, and the site may refuse the request; cache the blocks first
    with ``seamm-thermochemistry-atct covariance``).

The results that can be stored are the enthalpy of formation at 0 K and at
298.15 K, its 95 % reference uncertainty, the rung and the reaction used, the
ATcT version and, when the target is in ATcT, its values. They are also stored
as the properties ``DfH(0 K)``, ``DfH(298.15 K)`` and ``DfH uncertainty (95%)``
of the configuration.

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
properties are named with the Model Chemistry's level, for example
``H#Thermochemistry#MOPAC:SQM@PM6-ORG`` or
``S#Thermochemistry#ORCA:DFT@B3LYP/bse:def2-TZVP``, or with the subflowchart
engine's model, e.g. ``H#Thermochemistry#PM6-ORG``.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
