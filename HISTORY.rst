=======
History
=======
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
