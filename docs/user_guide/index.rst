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
species, with the 95 % uncertainty of the ATcT part (using the ATcT
covariances, fetched from atct.anl.gov when they are not cached). If the target
itself is in ATcT the difference is reported too.

Radicals and ions are handled: every fragment containing the radical or
charged atom keeps that spin or charge, and the ATcT entry is selected with
the same charge and multiplicity. The reactions follow from the SMILES of the
structure, so the system's charge and spin multiplicity must agree with it;
the step stops if they do not. Build a radical from its SMILES
(``[CH2]C`` for the ethyl radical), which sets both. A fragment takes the spin
multiplicity of its ground state: its unpaired electrons + 1, except the
singlets that a SMILES writes with two unpaired electrons (SiH\ :sub:`2`,
SiCl\ :sub:`2`, CF\ :sub:`2`, BH, ...) and the triplets it writes closed-shell
(O\ :sub:`2`, OCCO, cyclic C\ :sub:`3`, ...; seamm_thermochemistry's
``cbh.ground_state_multiplicity``). The target keeps its structure's
multiplicity, and the step says so when that is not the ground state.

The ATcT tables are not distributed with SEAMM: the seamm-thermochemistry
installer downloads the current version, and a calculation that finds none
installed downloads it then (see the seamm_thermochemistry documentation; a
newer version is added only on purpose, with ``seamm-thermochemistry-atct
download``). A rung is usable only when all its reference species are in
ATcT, and the step reports the rungs it could not use and why; the result
stored is the highest usable rung, and every rung goes into the table ``CBH``.
When no rung is usable but the target is itself in ATcT (methane, benzene, the
methyl radical, CO: every rung reproduces the target), its ATcT value is the
result, labelled "target is a reference species", and nothing is computed.
The spread of the enthalpies of formation over the rungs from CBH-1 up is
reported and stored as an indicator of the uncertainty: a large spread points
to trouble with the method or with the reference values.

Calculate
    ``enthalpy of formation`` (default) runs the CBH reactions as described here.
    ``the species only (for the fragment store)`` computes only the structure
    itself, exactly as a fragment of a reaction: optimized and its harmonic
    frequencies at the geometry level, then its energy at the energy level, with
    the properties and settings fingerprints that ``seamm-thermochemistry-fragments
    harvest`` files in the store. It runs no reactions and needs no ATcT entry, so
    it can add any species to the store, for example to recompute stored species
    after a change of the settings that the fingerprints record. The species keeps
    the structure's own name; harvest it with ``--species-map`` or ``--all-systems``
    (or name the systems ``CBH:<SMILES>``). Only the controls that apply are
    shown, and the ZPE is stored unscaled, as for every species.

Preset
    The levels of theory and rungs, chosen together:

    =========  ==========================  ==============================  ==========
    Preset     Geometries and frequencies  Energies                        Rungs
    =========  ==========================  ==============================  ==========
    quick      r2SCAN-3c                   r2SCAN-3c                       up to CBH-2
    default    r2SCAN-3c                   revDSD-PBEP86-D4/def2-TZVPP     up to CBH-2
    custom     the Model Chemistry         as set below                    as set below
    =========  ==========================  ==============================  ==========

    The presets compute every species with ORCA, so they need no Model Chemistry step
    before this one, and scale the ZPEs by 0.9899. Each uses the highest rung
    available, at most CBH-2: CBH-3 is no more accurate (its error is dominated by
    the reference values), and *Include CBH-3* adds it only to the spread between
    the rungs. One r2SCAN-3c geometry and frequency calculation serves every
    preset, so a species stored in the fragment store at one preset is not
    recomputed for another; only its energy is.

    On 46 molecules against ATcT or vetted reference values at CBH-2, ``default``
    has an RMS error of about 2.1 kJ/mol and ``quick`` about 2.2, but where CBH-2
    is not available ``quick`` is much weaker: several kJ/mol at CBH-1, and 7 kJ/mol
    MAE on small aromatics, whose rings are kept whole so that CBH-2 is trivial.
    ``quick`` therefore warns when the highest rung is below CBH-2, and suggests
    ``default``. (A ``high`` preset, CCSD(T)-F12/cc-pVTZ-F12 energies, about
    1 kJ/mol, is to follow.)

ZPE scale factor
    ``from the table`` (default): the factor for the geometry level, with its
    source, which the output and the results give (r2SCAN-3c 0.9899, fitted to
    small molecules; PBE0/def2-TZVP 0.975, a value for another basis;
    B2PLYP-D3/cc-pVTZ 0.9863). A level not in the table is not scaled, with a
    warning. Or a number.

CBH rungs
    ``all usable`` (default for ``custom``), or a list such as ``1, 2``.

Keep whole
    Groups kept as single units when the reactions are built, so that they are
    never broken: by default aromatic ring systems, 3- and 4-membered rings and
    charge-separated or hypervalent groups (nitro, sulfonyl, phosphoryl, PF6-,
    ...). ``none`` builds the plain hierarchy.

    Aromatic ring systems of any size and with any heteroatoms are kept whole,
    fused rings merged into one unit; non-aromatic 5- and 6-membered rings
    (cyclopentane, cyclohexane, THF, EC, FEC, ...) are opened like chains. This
    default is the most accurate: on 14 molecules whose reactions differ either
    way, the mean absolute error against reference values was 0.9 kJ/mol kept
    whole against 2.0 opened with F12 energies, 2.2-2.3 against 2.9-3.5 with
    revDSD, and 6.0 against 10.8 with PBE0-D4. Opening a 3-membered ring leaves
    its strain to the method (methylcyclopropane: 0.1 against 5.0 kJ/mol with
    F12), and opening benzene into polyene fragments at CBH-3 is 3 to 6 kJ/mol
    off in every method tried.

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

Energies at
    ``the geometry level`` (default): the Model Chemistry gives the energies as
    well as the geometries, frequencies, ZPE and thermal corrections. Or another
    model chemistry, for example ``ORCA:DFT@REVDSD-PBEP86-D4_2021/def2-TZVPP`` or
    ``ORCA:QC@CCSD(T)-F12/cc-pVTZ-F12``: its single point on each species'
    geometry replaces the energy, and each species' enthalpy becomes

        H(T) = E(energy level) + [H(T) - E](geometry level)

    so that a costly energy is used where it matters, on cheaper geometries and
    frequencies. The target's single point runs first; the fragments' then run
    as one batch of separate calculations, which can go to a cluster queue and
    can use methods that have no gradients (F12, for instance). This needs the
    computed species with harmonic thermal corrections. The species table shows
    both energies, and the results record both levels.

    With the fragment store in use, a fragment whose geometry and frequencies
    at the geometry level are in the store is not recomputed, and one whose
    energy at the energy level, with the same settings, is stored on that
    geometry is not run again. The energy-level energies, their settings
    fingerprints and the geometry each ran on are stored on the species, so
    that ``seamm-thermochemistry-fragments harvest`` files them in the store's
    geometry and energy tables.

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
    ``use and add what is computed`` (default), ``use only`` or ``ignore``.
    Species already in the fragment-energy store at this model chemistry are
    taken from it instead of being computed, and serve as the fallback in
    "existing results" mode. With ``use and add what is computed`` the species
    the step computes are added to the store at once, so later iterations of a
    Loop and later jobs find them; the output lists what was added and anything
    refused (a structure that rearranged while it was optimized, one with an
    imaginary mode). They are added only when the store is on a local disk of
    the machine running the job. Otherwise the job is marked, and is harvested
    into the store when it returns to its JobServer (or by hand with
    ``seamm-thermochemistry-fragments harvest``).

    The store is read in two layers: the installation's own (local) store
    first, then SEAMM's official one when it is installed. The species table
    says "store (official <version>)" for a value from the official store, and
    names both when an energy and its geometry come from different stores; the
    output notes any value on which the two disagree. Only the local store is
    ever written. *Fragment store file* names the local store. The species table says where
    each value came from; a record without a settings fingerprint is
    "unverified" and the output says so. When the engine provides them (ORCA),
    the settings fingerprints of each species' calculations are stored with its
    energies, so the harvest files verified records.

    Before the store, a fragment already computed in the same job -- by an
    earlier iteration of a Loop over many targets, say -- is taken from the
    job's own database ("job" in the species table), when it was computed with
    the same settings and has no imaginary mode; an energy-level single point
    is reused only on the geometry it ran on. Each species is looked up with
    the fingerprints expected for it, not the target's: ORCA uses exact
    exchange (NoCOSX) for atoms and for Na, Mg, Zn, B and P, so the methane
    fragment of a boron or phosphorus target has a different fingerprint from
    the target's. If a computed fingerprint is not the one predicted, the
    output says so; the species is still used.

ATcT version, directory, covariances
    The version of the Thermochemical Network (``latest`` cached by default),
    where the versions are kept, and whether to fetch missing covariance blocks
    from atct.anl.gov (``fetch if needed``, the default) or use only those
    cached (``cached only``). A block that is missing, or that cannot be fetched
    (a compute node often has no network, and the site may refuse the request;
    the step stops trying after the first failure), leaves the reference
    uncertainty in quadrature, without the correlations. The result's
    *uncertainty method* says which, and ``seamm-thermochemistry-atct
    covariance`` caches the blocks ahead of time.

The results that can be stored are the enthalpy of formation at 0 K and at
298.15 K, its 95 % reference uncertainty and how it was obtained, the rung and the reaction used, the
ATcT version and, when the target is in ATcT, its values. They are also stored
as the properties ``DfH(0 K)``, ``DfH(298.15 K)`` and ``DfH uncertainty (95%)``
of the configuration.

Every species the step computes, the target and each fragment, is written as
it was computed (optimized, if optimization is on) to ``<SMILES>.sdf`` in the
step's directory, and all of them together to ``species.sdf``, with the SMILES,
the role (target or fragment), the model chemistry, E, ZPE, H(0 K),
H(298.15 K) and the number of imaginary modes as SDF data fields. Open them to
check that each conformation is the one intended; the output also lists the
heavy-atom dihedrals of the final structure next to those of the starting one.
The files also hold the properties stored on the structure, and can be read
back with Read Structure, for example to restart a species from its computed
geometry. Slashes in a SMILES (double-bond stereo) become ``_`` and ``-`` in the
file name. Species taken from the fragment store or from existing results are not
written, since their structures are not computed here.

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
