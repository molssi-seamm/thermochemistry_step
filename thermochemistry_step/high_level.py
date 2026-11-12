# -*- coding: utf-8 -*-

"""The energy level of the 'high' preset: canonical CCSD(T)-F12 or DLPNO-CCSD(T)-F12,
chosen per target by thermo's rules (science, 2026-10-30, reaffirmed 2026-11-08 and
2026-11-11):

1. Canonical CCSD(T)-F12D/RI cc-pVTZ-F12 up to 546 basis functions (fits in 96 GB);
   above that it needs 200-330 GB.
2. Above 546: DLPNO-CCSD(T)-F12D with TightPNO and TCutPNO 1e-8, only for closed-shell,
   saturated C/H/O(/N) molecules, where CH2F/CHF2/CF3 on a non-fluorinated carbon is
   fine but two bonded fluorinated carbons are not (every DLPNO-F12 validation, jobs
   5516, 5536, 5560, was closed-shell).
3. pi systems (alkenes, dienes, aromatics, carbonyl/ester/carbonate) and perfluoro
   segments stay canonical at any size: above 546 they get no 'high' result.
4. DLPNO and canonical never mix within one reaction: the fragments follow the
   target's level (the step uses one energy level for the whole reaction).

Anything else above 546 basis functions gets no 'high' result, with the reason.
"""

#: The levels (orca_step's F12 model chemistries)
CANONICAL = "ORCA:QC@CCSD(T)-F12D_RI/cc-pVTZ-F12"
DLPNO = "ORCA:QC@DLPNO-CCSD(T)-F12D(TightPNO,TCutPNO=1e-8)/cc-pVTZ-F12"

#: The largest canonical CCSD(T)-F12 calculation (basis functions), rule 1
CANONICAL_MAX_BF = 546

#: ORCA 6.1.1's own cc-pVTZ-F12 basis functions per element (spherical), measured
#: from ORCA's output ("Number of basis functions") for each atom, 2026-11-11:
#: heptane is 7 x 60 + 16 x 18 = 708, as in job 5474. The Basis Set Exchange's
#: version is smaller (53 for C), so ORCA's counts are used. ORCA has cc-pVTZ-F12
#: for H-Ar only; other elements give no 'high' result.
NBF_CC_PVTZ_F12 = {
    1: 18, 2: 19,
    3: 60, 4: 60, 5: 60, 6: 60, 7: 60, 8: 60, 9: 60, 10: 60,
    11: 96, 12: 93,
    13: 90, 14: 90, 15: 90, 16: 90, 17: 90, 18: 90,
}  # fmt: skip

#: The elements DLPNO-F12 was validated for (rule 2)
DLPNO_ELEMENTS = {1, 6, 7, 8, 9}


def basis_functions(atomic_numbers):
    """ORCA's cc-pVTZ-F12 basis functions for the atoms, or (None, element) when an
    element's count is not known."""
    n = 0
    for Z in atomic_numbers:
        if Z not in NBF_CC_PVTZ_F12:
            return None, Z
        n += NBF_CC_PVTZ_F12[Z]
    return n, None


def dlpno_problem(smiles, multiplicity):
    """Why DLPNO-F12 is not validated for this molecule, or '' if it is."""
    if multiplicity != 1:
        return (
            f"it is open-shell (multiplicity {multiplicity}), and DLPNO-F12 is "
            "validated for closed-shell molecules only"
        )
    from seamm_thermochemistry.smiles import mol_from_smiles

    mol = mol_from_smiles(smiles)
    if mol is None:
        return f"its SMILES '{smiles}' cannot be read"
    from rdkit import Chem

    mol = Chem.AddHs(mol)
    others = sorted(
        {
            a.GetSymbol()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() not in DLPNO_ELEMENTS
        }
    )
    if others:
        return (
            f"it contains {', '.join(others)} (DLPNO-F12 is validated for C, H, O, "
            "N and F)"
        )
    for bond in mol.GetBonds():
        if bond.GetIsAromatic() or bond.GetBondType() != Chem.BondType.SINGLE:
            return (
                "it has a pi system (a multiple or aromatic bond), which needs "
                "canonical F12"
            )

    def fluorinated(atom):
        return atom.GetAtomicNum() == 6 and any(
            n.GetAtomicNum() == 9 for n in atom.GetNeighbors()
        )

    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if fluorinated(a) and fluorinated(b):
            return (
                "it has two bonded fluorinated carbons (a perfluoro segment), which "
                "need canonical F12"
            )
    return ""


def high_level(smiles, multiplicity, atomic_numbers, max_bf=CANONICAL_MAX_BF):
    """The 'high' energy level for a target, by the rules above.

    Parameters
    ----------
    smiles : str
        The target.
    multiplicity : int
    atomic_numbers : [int]
        All the atoms, hydrogens included.
    max_bf : int
        The largest canonical CCSD(T)-F12 calculation, basis functions (546 for
        96 GB nodes; more with more memory).

    Returns
    -------
    (str or None, str)
        The level (model chemistry), or None when there is no 'high' result for
        this target; and the reason, for the output.
    """
    n, unknown = basis_functions(atomic_numbers)
    if n is None:
        return None, (
            f"ORCA has no cc-pVTZ-F12 basis for element {unknown} (it has H-Ar)"
        )
    if n <= max_bf:
        return CANONICAL, (
            f"{n} cc-pVTZ-F12 basis functions (at most {max_bf}): "
            "canonical CCSD(T)-F12"
        )
    problem = dlpno_problem(smiles, multiplicity)
    if problem:
        return None, (
            f"{n} cc-pVTZ-F12 basis functions is too many for canonical CCSD(T)-F12 "
            f"(at most {max_bf}), and DLPNO-CCSD(T)-F12 is not used because "
            f"{problem}"
        )
    return DLPNO, (
        f"{n} cc-pVTZ-F12 basis functions (over {max_bf}), a closed-shell "
        "saturated molecule: DLPNO-CCSD(T)-F12 (TightPNO, TCutPNO 1e-8)"
    )
