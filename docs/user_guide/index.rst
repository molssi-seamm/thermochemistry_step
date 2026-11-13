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
structure, so the system's charge must agree with it, and so must its number of
electrons: the step stops if they do not. Build a radical from its SMILES
(``[CH2]C`` for the ethyl radical), which sets both. The structure's spin
multiplicity may be lower than the unpaired electrons the SMILES shows: SiCl\ :sub:`2`
and BH are singlets that a SMILES writes with two unpaired electrons. The step
then uses the structure's multiplicity and says so, noting whether it is the
known ground state. A fragment takes the spin
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
The target is found in ATcT also when ATcT writes it as another resonance form:
SiO from a structure is ``[O+]#[Si-]`` but ATcT's is ``[Si]=O``, and OBrO is
``[O][Br+][O-]`` against ``[O]Br[O]``. The match ignores bond orders, formal
charges and radical centres but keeps the atoms, their connections and
hydrogens, and the total charge, so isomers such as HCN and HNC stay apart; it
is used only when exactly one ATcT species matches, and the output names the
ATcT SMILES it matched.
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
    high       r2SCAN-3c                   CCSD(T)-F12/cc-pVTZ-F12         up to CBH-2
    custom     the Model Chemistry         as set below                    as set below
    =========  ==========================  ==============================  ==========

    The presets compute every species with ORCA, so they need no Model Chemistry step
    before this one, and scale the ZPEs by 0.9899. ``quick`` and ``default`` are the
    composite recipes ``SEAMM:CMP@r2SCAN-3c-HO`` and
    ``SEAMM:CMP@revDSD-r2SCAN-3c-HO`` when the Composite step (composite-step) is
    installed, and the results record the recipe, its version and its fingerprint;
    the levels and the numbers are the same either way. A composite recipe may also
    be the flowchart's model chemistry with ``custom``, for the CBH approach only,
    or the energy level. A recipe that combines several calculations runs each
    component as its own batch on every species; each component is looked up on its
    own in the job and in the fragment store, and stored as its own energy, never the
    combination. The species table names how each component was obtained, and a
    species whose component fails stops the step with that component named.

    Every energy the step computes at the energy level, a component or not, is
    stored with its parts (HF, correlation, triples, ...) when the program reports
    them. A recipe that uses parts takes them from the job or the store with the
    energy; a stored energy without parts is computed again, and the store then
    keeps its parts.

    Each preset uses the highest rung
    available, at most CBH-2: CBH-3 is no more accurate (its error is dominated by
    the reference values), and *Include CBH-3* adds it only to the spread between
    the rungs. One r2SCAN-3c geometry and frequency calculation serves every
    preset, so a species stored in the fragment store at one preset is not
    recomputed for another; only its energy is.

    When the energies are at another level than the geometries (``default``,
    ``high``, or ``custom`` with an energy level), the step also gives the same
    reactions at the geometry level alone: each species' r2SCAN-3c energy with the
    same thermal corrections and reference values. That is the ``quick`` result,
    so a ``default`` run need not be repeated with ``quick``. The output prints its
    rung table after the main one, and the results ``DfH 0K (geometry level)``,
    ``DfH 298K (geometry level)``, ``DfH uncertainty (geometry level)`` and
    ``CBH rung spread (geometry level)`` hold it, with ``CBH geometry level
    preset`` naming the preset it equals (``quick``, or empty for other geometry
    levels).

    On 46 molecules against ATcT or vetted reference values at CBH-2, ``default``
    has an RMS error of about 2.1 kJ/mol and ``quick`` about 2.2, but where CBH-2
    is not available ``quick`` is much weaker: several kJ/mol at CBH-1, and 7 kJ/mol
    MAE on small aromatics, whose rings are kept whole so that CBH-2 is trivial.
    ``quick`` therefore warns when the highest rung is below CBH-2, and suggests
    ``default``.

    ``high`` takes the energies from CCSD(T)-F12, about 1 kJ/mol, chosen for each
    target:

    * canonical CCSD(T)-F12D/RI with cc-pVTZ-F12 when the target has at most 546
      basis functions (*Largest canonical F12 (basis functions)*; 546 fits in
      96 GB of memory, and larger calculations need 200-330 GB);
    * otherwise DLPNO-CCSD(T)-F12D (TightPNO, TCutPNO 1e-8), but only for
      closed-shell molecules of H, C, N, O and F with no multiple or aromatic bonds
      and no two bonded fluorinated carbons, as validated;
    * otherwise no ``high`` result, with the reason, and the ``default`` result
      instead.

    The fragments follow the target's choice, so canonical and DLPNO energies are
    never mixed in one reaction. ORCA has cc-pVTZ-F12 for H to Ar only. ``high``
    needs orca-step 2026.11.11.1 or later.

ZPE scale factor
    ``from the table`` (default): the factor for the geometry level, with its
    source, which the output and the results give (r2SCAN-3c 0.9899, fitted to
    small molecules; PBE0/def2-TZVP 0.975, a value for another basis;
    B2PLYP-D3/cc-pVTZ 0.9863). A level not in the table is not scaled, with a
    warning. Or a number. A level with a frequency correction (below) is not
    scaled: the two are never combined, and a number given with a correction is
    refused.

Frequency correction
    ``from the table`` (default) or ``none``. Some methods get their frequencies
    wrong in a way one scale factor cannot fix, and have a correction applied to
    each frequency before the ZPE, H(T) - H(0), S and Cp are computed, in the
    harmonic approach and for the CBH species alike. So far there is one, for PM7
    (Rozanska et al., J. Chem. Eng. Data 59, 3136 (2014), eq 1): unchanged below
    125 cm\ :sup:`-1`, 1.006 ν - 27.7 up to 2500 cm\ :sup:`-1`, and 3.99 ν - 7849.2
    above. Other MOPAC methods (PM6, PM6-ORG, ...) have none: PM7's fit does not
    carry over, and the output says so.

    The output lists the harmonic and corrected frequencies, the ZPE with and
    without the correction, and warnings:

    * frequencies near the switches, where the map is discontinuous (raw
      98-125 and 2400-2600 cm\ :sup:`-1` for PM7: 124 stays 124 but 126 becomes
      99.1; 2499 becomes 2486.3 but 2501 becomes 2129.8);
    * stretches of hydrogen bonded to elements outside the fit set (C, H, N and O
      for PM7), found from the normal modes;
    * frequencies above the range of the fit, which are left uncorrected (PM7:
      above a raw 2900 cm\ :sup:`-1`, where experiment ends at about 3700; its H2
      stretch would otherwise become 14400 cm\ :sup:`-1`). The species is marked
      ``partly not applicable``.

    Without thermal corrections (``none (the energy is the enthalpy)``) the
    frequencies are not used and the correction does nothing, which the output
    says. The correction's tag goes into the frequency fingerprint of the
    species' thermal data, so corrected and uncorrected values are never mixed in
    the fragment store.

Correct
    With a frequency correction: ``full`` (default), the correction as published,
    or ``high-frequency only``, leaving out the branches the source calls
    optional (for PM7, the 0.6 % change below 2500 cm\ :sup:`-1`).

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
    imaginary mode), and names each computed species that was not added, with
    the reason. Without thermal corrections (MOPAC's enthalpy of formation as
    the enthalpy) only the energy is stored, with no ZPE or enthalpy, so a run
    with harmonic thermal corrections at the same level computes its own. They
    are added only when the store is on a local disk of
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

    The target, too, is taken from the store (*Reuse a stored target*, ``if the
    structure matches`` by default) when the store has it at this geometry level
    and the stored structure is the one the target starts from: the same atoms in
    the same order, the same bonds, and the heavy atoms within 0.3 Å RMSD after
    the best superposition, with no mirroring, so another conformer or stereoisomer
    is never taken. The output gives the RMSD, and the target's coordinates become
    the stored ones; its energy at another energy level is then looked up in the
    store or computed on that geometry. A stored target that does not match is
    reported, and the target is computed. ``never`` always computes it.

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

Supplementary references, Use my own references
    Reference species that ATcT lacks (whole elements such as Li, Na, P or Al,
    and many organic parents) can come from the reference library, which is
    installed beside the ATcT data. It holds evaluated compilations and
    measurements: printed JANAF, CODATA, the NBS tables of 1982, spectroscopic
    D\ :sub:`0`, W4-17 and measurements from the NIST WebBook, each cited to its
    original source, and computed references.

    - ``compilations and computed``, the default, uses all of them;
      ``compilations`` leaves out the computed references; ``none`` uses ATcT
      only.
    - ATcT always comes first. The library is used only with ATcT as the
      reference values, not with a ladder.
    - For an element that ATcT lacks, all its references are taken from one
      source where possible, so that the source's own error cancels in the
      reaction. When the references of one element come from several sources,
      the output says so and the result is labelled *mixed: several sources*.
    - A library value has no covariance with ATcT. Its uncertainty (95 %, plus
      the source's typical offset from ATcT) is added in quadrature.
    - A value tabulated at only one temperature is converted to the other with
      the species' own H(298.15 K) − H(0) from this calculation. The species
      table marks the value with ``*``.
    - *Use my own references* (shown unless the setting is ``none``) also uses
      the references you added to your local library with
      ``seamm-thermochemistry-references``, ahead of the compilations.

    A target that ATcT lacks but the library has is treated like a target in
    ATcT. When a rung works, the CBH value is the result, and its difference
    from the library value is reported (*CBH - library*) and flagged when it
    exceeds twice the library's uncertainty: a check of the library value. When
    no rung works, the library value is the result, unchecked.

    The result's *CBH references* says where the references of the reaction
    came from: ``ATcT-only``, ``mixed``, ``mixed: several sources``, ``computed
    references`` or ``user references``. *references* records the versions and
    the library values used.

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

Each computed species is checked after its optimization:

* An optimization that has not converged in *Fragment max steps* continues from
  where it stopped for as many steps again. One that still has not converged is
  used, but listed in the result ``CBH unconverged species`` (empty when all
  converged), marked ``optimization converged: no`` in its SDF file, and not added
  to the fragment store.
* A target whose bonds are no longer those of its SMILES after the optimization
  (dimethyl sulfate falling apart to SO\ :sub:`2` and CH\ :sub:`3`\ OOCH\ :sub:`3`,
  say) fails, since its energy would be that of other molecules. A fragment that
  does so is used, listed in ``CBH changed fragments`` and not added to the store.

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
