=======
History
=======
2026.11.12.2 -- Reference species that ATcT lacks, from the reference library
    * New: the CBH approach can take reference species that ATcT lacks from the
      supplementary reference library of seamm-thermochemistry (printed JANAF,
      CODATA, the NBS tables of 1982, spectroscopic D0, W4-17, NIST WebBook
      measurements, computed references), with one source per element where
      possible. Their uncertainties (95 %, with each source's offset from ATcT) are added in
      quadrature, and a value given at one temperature only is converted with the
      species' own computed H(298.15 K) - H(0).
    * New parameters: "Supplementary references" (compilations and computed /
      compilations / none) and "Use my own references". The default,
      "compilations and computed", is Paul's decision 4 (2026-10-29). Flowcharts
      saved before this release take the default when they are run again, so a
      target that had no usable rung with ATcT alone may now get a result from
      mixed references. Targets that ATcT alone closes are unchanged, and "none"
      reproduces the earlier results exactly.
    * New results: "CBH references" (ATcT-only, mixed, mixed: several sources,
      computed references, user references) and "references" (the versions of
      ATcT and the libraries, and the library values used).
    * A target that ATcT lacks but the library has: with a usable rung, the CBH
      value is checked against the library value ("CBH - library", flagged beyond
      twice the library's uncertainty); with none, the library value is the
      result.
    * Needs seamm-thermochemistry 2026.11.12.1 or later.

2026.11.12.1 -- The parts of energies kept and reused
    * Every energy the CBH approach computes at the energy level is stored with its
      parts (HF, correlation, triples, ...) when the program reports them, so
      composite recipes that use the parts find them in the job or the fragment
      store.
    * A stored energy without the parts a recipe needs is computed again, and the
      store then keeps its parts.
    * Needs seamm-thermochemistry 2026.11.12 or later.

2026.11.12 -- Composite recipes with several calculations as the CBH energy level
    * A composite recipe that combines several calculations (``SEAMM:CMP@<recipe>``)
      can be the energy level of the CBH approach, as the flowchart's model
      chemistry or as the energy model chemistry. Each component runs as its own
      batch on every species.
    * Each component is looked up on its own in the job and in the fragment store,
      and stored there as its own energy under its own fingerprint, never the
      combination, so later recipes that share a component reuse it. The
      composite's record goes with the results.
    * The species table names how each component was obtained. A species whose
      component fails stops the step, with the component named.
    * A rerun that takes every component from the store reproduces the composite
      energies to about 1e-12 E_h (not bit for bit: the store keeps energies through
      a unit conversion).
    * Needs composite-step 2026.11.12 and seamm-exec 2026.11.12 or later for recipes
      with several components.

2026.11.11.6 -- The 'high' preset (CCSD(T)-F12); presets as composite recipes
    * The new ``high`` preset of the CBH approach takes the energies from CCSD(T)-F12
      with cc-pVTZ-F12 on r2SCAN-3c geometries: canonical CCSD(T)-F12D/RI up to 546
      basis functions (a parameter; 546 fits in 96 GB), and above that
      DLPNO-CCSD(T)-F12D (TightPNO, TCutPNO 1e-8) for the closed-shell, saturated
      H/C/N/O/F molecules it was validated on. Other targets get no ``high`` result,
      with the reason, and the ``default`` preset's energies instead. The fragments
      follow the target, so canonical and DLPNO energies are never mixed in one
      reaction. Needs orca-step 2026.11.11.1 or later.
    * With the Composite step (composite-step) installed, the ``quick`` and
      ``default`` presets are the composite recipes ``SEAMM:CMP@r2SCAN-3c-HO`` and
      ``SEAMM:CMP@revDSD-r2SCAN-3c-HO``, and the results record the recipe, its
      version and its fingerprint. The levels and the numbers are unchanged. A
      composite recipe may also be the flowchart's model chemistry for the CBH
      approach.

2026.11.11.5 -- Bugfix: CBH targets in another resonance form, and singlet diradicals
    * A target that ATcT writes as another resonance form is now found: SiO from a
      structure is [O+]#[Si-] but ATcT's is [Si]=O, and OBrO is [O][Br+][O-] against
      [O]Br[O]. Such a target got no result when every rung reproduced it. The match
      ignores bond orders, formal charges and radical centres but keeps the atoms,
      their connections and hydrogens, and the total charge, so isomers such as HCN
      and HNC stay apart. It is used only when exactly one ATcT species matches, and
      the output names the ATcT SMILES it matched.
    * A structure's spin multiplicity may now be lower than the unpaired electrons its
      SMILES shows: SiCl2 and BH are singlets that a SMILES writes as diradicals, and
      they stopped the step. The structure's multiplicity is used and the output says
      so, noting whether it is the known ground state. The charge and the number of
      electrons must still agree with the SMILES.

2026.11.11.4 -- ATcT is downloaded when needed; covariances are fetched by default
    * A calculation that finds no ATcT version installed downloads the current one
      from atct.anl.gov (an installation made before the seamm-thermochemistry
      installer set ATcT up), and says so; parallel jobs may do so at once. If the
      download fails the error names the installer and the offline route. A
      specific version that is not installed is never downloaded.
    * The ATcT covariances default to "fetch if needed": a block that is not cached
      is fetched and cached. After the first failed fetch (e.g. no network) the
      step stops trying for the rest of the calculation, and those reactions'
      uncertainties are in quadrature, labelled "quadrature (covariance fetch
      failed)". Flowcharts saved earlier keep the setting they were saved with.
    * The uncertainty method is a result too ("DfH uncertainty method"), as well as
      a column of the CBH table.
    * Needs seamm-thermochemistry 2026.11.11.11 or later.

2026.11.11.3 -- CBH: reference species as targets; ground-state multiplicities
    * A target that is itself an ATcT species (methane, benzene, pyridine, the methyl
      radical, CO, ...), for which no rung gives a usable reaction, returns its ATcT
      value, labelled "target is a reference species", instead of stopping with "No
      CBH rung has all its reference species in ATcT". Nothing is computed.
    * The fragments, and their ATcT lookups, take the spin multiplicity of their
      ground state: the singlet SiH2, SiCl2, CF2, CCl2, BH, ... that a SMILES writes
      as triplets, and the triplet O2, S2, SO, OCCO and cyclic C3 that it writes
      closed-shell. The target keeps its structure's multiplicity, and the step says
      so when that is not the ground state.
    * Needs seamm-thermochemistry 2026.11.11.7 or later.

2026.11.11.2 -- The official fragment store is read too
    * The species are looked up in the installation's own (local) fragment store
      first, then in SEAMM's official store when it is installed. The species table
      says which store a value came from, and names both when an energy and its
      geometry come from different stores; disagreements between the two are
      noted in the output. The lowest conformer of either store is used, a tie
      going to the local one. Only the local store is ever written.
    * Needs seamm-thermochemistry 2026.11.11.2.

2026.11.11.1 -- The species computed are added to the fragment store
    * "Fragment store" has the choices "use and add what is computed" (the new
      default), "use only" (formerly "use if available") and "ignore". With the
      default, each species the isodesmic (CBH) approach computes -- the target, the
      fragments, or the species alone -- is added to the fragment store as soon as
      it is done, with its geometry, frequencies and energies, so that later
      iterations of a Loop and later jobs reuse it. The output lists what was added
      and what was refused: a structure whose bonds are no longer those of its
      species, or one with an imaginary mode.
    * The species are added only when the store is on a local disk of the machine
      running the job. Otherwise the job is marked, and the JobServer harvests it
      into the store when the job comes back to it.
    * Needs seamm-thermochemistry 2026.11.11.1.

2026.11.11 -- Fragments reused within a job and looked up by their own fingerprints
    * A fragment already computed in the same job, e.g. by an earlier iteration of
      a Loop over many targets, is now taken from the job's database instead of
      being optimized and its Hessian computed again, when its settings
      fingerprints are the expected ones and it has no imaginary mode. Its energy
      at the energy level is reused too, if it ran on the same geometry. A serial
      Loop over 177 compounds (SEAMM job 5684) had computed some 500 fragments a
      second time or more (H2 114 times, CH4 74).
    * Each species is looked up, in the job and in the fragment store, with the
      fingerprints its own calculation would have. The target's were used for
      all, but ORCA uses exact exchange (NoCOSX) for atoms and for Na, Mg, Zn, B
      and P: the ordinary fragments of a boron or phosphorus target were never
      found in the store, and were computed again in every such job.
    * Bugfix: after computing a fragment the step looked the following ones up
      with that fragment's fingerprints rather than the target's, so after a
      boron fragment every later lookup in the iteration missed.
    * If a computed fingerprint is not the one predicted for the species, the
      output says so. The species is still used.
    * Needs seamm-thermochemistry, seamm-exec and orca-step 2026.11.11.

2026.11.10 -- The species only, for the fragment store
    * The isodesmic (CBH) approach can compute only the structure itself
      ("Calculate: the species only (for the fragment store)"), exactly as a
      fragment of a reaction: optimized with its Hessian at the geometry level,
      its energy at the energy level, and the properties and settings
      fingerprints that the fragment store's harvest reads. It runs no reactions
      and needs no ATcT entry, so any species can be added to the store, e.g. to
      recompute stored species at the settings the presets look up.

2026.11.8.1 -- Presets for the isodesmic approach; scaled ZPEs; the rung spread
    * The isodesmic (CBH) approach has presets: 'quick' (r2SCAN-3c), and 'default'
      (r2SCAN-3c geometries and frequencies, revDSD-PBEP86-D4/def2-TZVPP energies), the
      new default, both with ZPEs scaled by 0.9899 and the highest rung available up to
      CBH-2 (CBH-3 on request). They need no Model Chemistry step. 'custom' keeps the
      earlier behaviour. A flowchart that uses the CBH approach and was saved before
      this version gets the 'default' preset; choose 'custom' to keep its Model
      Chemistry.
    * The ZPEs are scaled by the factor for the geometry level from a table with its
      source (r2SCAN-3c, PBE0/def2-TZVP, B2PLYP-D3/cc-pVTZ), or by a given factor; the
      factor and its source are printed and stored.
    * The spread of the enthalpy of formation over the rungs from CBH-1 up is reported
      and stored as an indicator of the uncertainty, and 'quick' warns when the
      highest rung is below CBH-2.

2026.11.8 -- Energies at another level than the geometries
    * The isodesmic (CBH) approach can take the species' energies from a second model
      chemistry ("Energies at"), e.g. revDSD-PBEP86-D4/def2-TZVPP or
      CCSD(T)-F12/cc-pVTZ-F12, on the geometries and with the frequencies of the
      Model Chemistry: H(T) = E(energy level) + [H(T) - E](geometry level). The
      single points run as separate calculations, the fragments' as one batch, so
      they can go to a cluster queue and can use methods without gradients.
    * With the fragment store, stored geometries and energies are reused, and the
      new energies, frequencies, symmetry numbers and moments of inertia are stored
      on the species for the store's harvest.

2026.11.5 -- Internal: the Hessian driver moved to seamm_exec; the CBH ring rules documented
    * The code that drives the Model Chemistry for energies, optimizations and Hessians
      moved to seamm_exec (2026.11.5), where the Normal Mode Sampling and Vibrational
      Analysis steps share it. Nothing changes in how this step runs.
    * The user guide now says which rings the CBH reactions keep whole (aromatic ring
      systems, rings of up to 4 atoms, charge-separated groups) and which they open
      (non-aromatic 5- and 6-membered rings), with the evidence for that default.

2026.10.30 -- Bugfix: the species SDF files can be read back
    * The ``<SMILES>.sdf`` and ``species.sdf`` files of the Isodesmic (CBH) approach
      lacked the blank line that ends the last of the structure's stored properties, so
      a reader joined that value with the next data field and Read Structure could not
      read the file ("could not convert string to float"). The files are now valid SDF
      and can be read back, for example to restart a species from its computed
      geometry. Files written by 2026.10.29.2 can be repaired by adding a blank line
      before ``> <SMILES>``.
    * In the output, each species' header now comes before the line saying which
      starting conformer was used; before, the target's conformer line was printed
      above its header.

2026.10.29.2 -- The structure of every computed CBH species is written to SDF files
    * Each species the Isodesmic (CBH) approach computes, the target and every fragment,
      is written as computed (optimized) to ``<SMILES>.sdf`` in the step's directory,
      and all of them to ``species.sdf``, with the SMILES, the role, the model
      chemistry, E, ZPE, H(0 K), H(298.15 K) and the number of imaginary modes as SDF
      data fields, so the conformations can be checked.
    * The output gives the heavy-atom dihedrals of each final structure next to those
      of the starting conformer.

2026.10.29.1 -- Bugfix: the MDI engine of a module-installed code now gets its modules
    * A program installed with environment modules (``installation = modules`` in its
      ``.ini`` file, as ORCA is on TinkerCliffs) is now started as a warm MDI engine with
      the environment that ``module load`` gives it. Before, the engine was started
      without the modules, so ORCA on TinkerCliffs could not load its own libraries and
      every Isodesmic (CBH) species failed at its first energy with "MDI_Recv failed".
      A normal (non-MDI) run of the program was not affected.

2026.10.29 -- Bugfix: MOPAC (and other single-threaded engines) could not be started
    * Since 2026.10.28.7 the step passed the core count to every program's MDI engine
      command, which only ORCA's accepts; MOPAC's engine then failed to start ("got
      an unexpected keyword argument 'ncores'"). The core count is now passed only to
      a program whose engine command takes it.

2026.10.28.10 -- Bugfix: one engine rank per physical core
    * Where hyperthreading is on, Slurm counts hardware threads as CPUs and confines
      the job to them, so a 16-task job holds 8 physical cores and the engine's 16
      ranks ran two to a core. The engine now uses one rank per distinct physical
      core among the CPUs the job may use (ask Slurm for twice the tasks, or
      ``--hint=nomultithread``). The pin to Slurm's abstract CPU ids of 2026.10.28.9,
      which could not apply, is gone.

2026.10.28.9 -- Bugfix: with the engine route the Hessian still went to tasks
    * When the Model Chemistry was driven as the warm engine (the CBH approach's
      default) and the program prefers separate tasks (ORCA), the structure was
      optimized over the engine but the finite-difference Hessian went to tasks, and
      no settings fingerprints were stored. The engine route now computes the finite
      differences over the engine as well, and the fingerprints are kept whichever
      branch runs. Under Slurm the engine and its MPI ranks are pinned to the job's
      allocated CPUs, so that two ranks do not share a physical core on a node
      without cgroups.

2026.10.28.8 -- Bugfix: parallel ORCA as the MDI engine could not find OpenMPI
    * The engine's own ORCA runs did not get the OpenMPI runtime the ORCA step gives
      its runs (``library-path`` in orca.ini, its ``mpirun`` first on the PATH, the
      shared-memory transport), so a parallel engine died in start-up ("cannot open
      libmpi.so.40"). The step now sets that environment, from the program's step,
      while the engine lives, with MPI core binding off, and restores it afterwards.

2026.10.28.7 -- Bugfix: the MDI engine ran on one core
    * The warm MDI engine (ORCA) was launched with one core whatever the job had. It
      now gets the job's allocation under Slurm (else the CPUs this process may use),
      capped by the job's ``ncores`` option when that is a number, and the step says
      how many. The engine's MPI runs are told not to pin their ranks to cores
      (``OMPI_MCA_hwloc_base_binding_policy=none``), as the ORCA step does.

2026.10.28.6 -- CBH species start from their lowest conformer
    * Each fragment built from its SMILES now starts from the lowest MMFF94 conformer
      of several RDKit embeddings ("Fragments start from", 20 embeddings by default),
      and the target can be rebuilt the same way ("Target starts from"), as CBH-ANL
      starts every species from its lowest rotamer: a single embedding of pentane to
      octane is gauche, and an optimization stays in that basin. The chosen chain
      dihedrals are printed. Requires seamm-thermochemistry 2026.10.28.4.

2026.10.28.5 -- Bugfix: the CBH species were not optimized when the model prefers tasks
    * For a model chemistry with an analytic gradient but no analytic Hessian whose
      program prefers separate tasks (ORCA with a double hybrid, for one) the step
      took the tasks route, so the CBH species were never optimized and no settings
      fingerprints were stored. A new parameter, "Energies and Hessian via", chooses
      the route: automatic (as before), the warm MDI engine, or separate tasks. The
      CBH approach now drives the engine by default, so every species is optimized
      over it and the Hessian is finite differences over the same connection; it
      falls back to tasks, with a note, only when the program has no engine or is
      not installed where the job runs.

2026.10.28.4 -- Settings fingerprints of the CBH species from the engine
    * After each CBH species the step asks the MDI engine for the settings
      fingerprints of its calculations (seamm-mdi 2026.10.28; ORCA's engine first) and
      stores them with the species' energies as the properties ``settings
      fingerprint``, ``element fingerprints`` and ``calculation``, and the ``frequency
      ...`` three for the Hessian (the analytic-Hessian record, or the energy record for
      finite differences), so that ``seamm-thermochemistry-fragments harvest`` files
      them as verified records and the fragment store is searched by fingerprint.
      Engines that give no fingerprints leave the records unverified, as before.

2026.10.28.3 -- CBH: radicals and ions, the fragment store, laddered reference values
    * Radicals and ions in the CBH approach: every fragment keeps the spin or charge of
      its atoms, the ATcT entry is selected with the same charge and multiplicity, and
      each species runs through the engine with its own charge and multiplicity. The
      target's SMILES is checked against the system's charge and spin multiplicity
      before the reactions are built (needs molsystem 2026.10.28 for radicals).
    * "Fragment store": species already in the fragment-energy store at this model
      chemistry (``seamm-thermochemistry-fragments harvest`` promotes them from finished
      jobs) are used instead of being computed, and serve as the fallback in "existing
      results" mode; the species table says where each value came from.
    * "Reference values": ATcT for every reference species, a laddered reference set
      built with ``seamm-thermochemistry-ladder`` (CBH-ANL), or ATcT below an
      uncertainty threshold and the ladder otherwise. With a ladder the uncertainty is
      the ATcT part through the ladder's coefficients on the anchors plus the ladder's
      measured method term, both reported.
    * Bugfix: the "existing results" mode could not read the stored energies (the
      property values came back as dicts).
    * Requires seamm-thermochemistry 2026.10.28.3.

2026.10.28.2 -- Enthalpies of formation by the isodesmic (CBH) approach with ATcT
    * A second approach, "Isodesmic (CBH) with ATcT references", gives the enthalpy
      of formation of a molecule at 0 K and 298.15 K from the connectivity-based
      hierarchy (CBH) of error-cancelling reactions, with the reference species taken
      from the Active Thermochemical Tables (ATcT) and the 95 % uncertainty of the
      ATcT part, using the ATcT covariances when cached. The target and every
      reference species are optimized with the engine and their harmonic
      thermochemistry calculated over the same warm connection; alternatively the
      energies are read from properties already stored on the systems
      ``CBH:<SMILES>``. For MOPAC and other NDDO methods, whose energy is already an
      enthalpy of formation, "Thermal corrections: none" skips the Hessians.
    * Every usable rung goes into the table ``CBH``; the stored results are those of
      the highest usable rung, with the difference from ATcT when the target is in
      the tables. The dialog shows only the settings that apply to the chosen
      approach and mode.
    * Requires seamm-thermochemistry 2026.10.28.2 and RDKit, and a cached ATcT
      version (``seamm-thermochemistry-atct download``); ATcT data are never
      distributed with SEAMM.

2026.10.28.1 -- The energies and Hessian come from the Model Chemistry
    * The step now takes its energies and gradients from the flowchart's Model
      Chemistry by default, driving the program as a warm engine over MDI instead of
      running a subflowchart once per displaced structure (a program start-up each
      time). When the program provides an analytic Hessian (ORCA for HF and most
      DFT functionals) it is used directly; otherwise the Hessian is finite
      differences of the gradients over the same connection, or as separate tasks
      where the job runs its tasks. The route is discovered at run time; there is
      nothing to set.
    * The subflowchart remains available ("Energies and gradients from:
      subflowchart") for calculations an MDI engine cannot do, and a flowchart from
      before this change, with a subflowchart and no Model Chemistry, still runs as
      it did. The dialog shows the Flowchart tab and the subdirectory handling only
      for the subflowchart, and reminds you to add a Model Chemistry step otherwise.
    * The properties are named with the Model Chemistry's level, e.g.
      ``H#Thermochemistry#ORCA:QC@HF/bse:def2-SVP``.
    * Checked: water with MOPAC PM6-ORG gives the same thermochemistry both ways
      (ZPE 49.76 vs 49.75 kJ/mol, H-E 59.680 vs 59.675 kJ/mol at 298.15 K), and
      ORCA HF/def2-SVP gives the same from its analytic Hessian and from finite
      differences (ZPE 60.67 kJ/mol both ways).

2026.10.28 -- Bugfix: default parameters crashed, ZPE counted twice, results not stored
    * The default symmetry number, "from system", stopped the step with
      ``NotImplementedError``. The rotational symmetry number is now perceived from the
      structure by testing rotation axes against the geometry (C60 gives 60, SF6 24,
      methane 12, benzene 12, water 2, a linear centrosymmetric molecule 2), and the
      step reports the number and the axes it found. Optimized geometries are symmetric
      to well within the 0.1 Å tolerance used.
    * The enthalpy and Gibbs energy corrections added the ZPE a second time, and in
      mixed units: ASE's H and G already include the ZPE, and the step added the ZPE in
      kJ/mol to values in kcal/mol. The ZPE also included the nominally zero
      translations and rotations. The thermal corrections U-E, H-E and G-E (each
      including the ZPE) and S are now in kJ/mol and J/mol/K, as the metadata says,
      and the ZPE comes from the true vibrations only. Water with PM6-ORG gives
      H-E-ZPE = 9.92 kJ/mol at 298.15 K, i.e. 4RT, and S = 188.2 J/mol/K.
    * The results advertised as total energy, H(0 K), U, H, S and G were never set, so
      they could be picked in the results tab but nothing was stored. The step now keeps
      the energy of the undisplaced structure from the subflowchart and stores all of
      them; U, H, S and G are lists over the pressures and temperatures, like the
      thermal corrections, and are named with the engine's model chemistry, e.g.
      ``H#Thermochemistry#PM6-ORG``.
    * A single atom and a single temperature or pressure crashed; both work now.
    * The number of frequencies to print was compared against the wrong list.
    * CI now installs the package with uv from its declared dependencies; the conda
      test environment file is gone.

2026.3.1 -- Internal: switching from deprecated library pkg_resources to importlib

2025.1.31: Add properties and other results.
    * Added ZPE, H, G, S, P, T, number of saddle point modes, and the transition states
      frequency (for transition states) to the results.
      
2024.12.22: First real release.

2024.10.17 (2024-10-17)
-----------------------

* Plug-in created using the SEAMM plug-in cookiecutter.
