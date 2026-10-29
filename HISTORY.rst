=======
History
=======
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
