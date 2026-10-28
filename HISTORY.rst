=======
History
=======
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
