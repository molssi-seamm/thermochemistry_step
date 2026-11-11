# -*- coding: utf-8 -*-

"""The 'high' preset's choice of CCSD(T)-F12 level (thermo's rules 1-4)."""

import pytest

pytest.importorskip("rdkit")
pytest.importorskip("seamm_thermochemistry")

from rdkit import Chem  # noqa: E402

from thermochemistry_step import high_level as H  # noqa: E402


def choose(smiles, multiplicity=1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return H.high_level(
        smiles, multiplicity, [a.GetAtomicNum() for a in mol.GetAtoms()]
    )


def test_orcas_counts_reproduce_heptane():
    """ORCA 6.1.1: CH4 132 and H2 36 cc-pVTZ-F12 functions; heptane 708 (job 5474)."""
    assert H.basis_functions([6, 1, 1, 1, 1])[0] == 132
    assert H.basis_functions([1, 1])[0] == 36
    assert H.basis_functions([6] * 7 + [1] * 16)[0] == 708


@pytest.mark.parametrize(
    "smiles", ["CCCC", "c1ccccc1", "CC(=O)OC", "FC(F)(F)C(F)(F)F", "CCCCC"]
)
def test_small_targets_are_canonical(smiles):
    level, reason = choose(smiles)
    assert level == H.CANONICAL and "canonical" in reason


def test_large_saturated_closed_shell_is_dlpno():
    # pentane 498 is canonical; hexane 612 is over 546
    assert choose("CCCCCC")[0] == H.DLPNO
    level, reason = choose("CCCCCCCC")  # octane, 8 x 60 + 18 x 18 = 804
    assert level == H.DLPNO and "804" in reason
    assert choose("OCCCCCCCO")[0] == H.DLPNO
    # CF3 on a non-fluorinated carbon is fine (rule 3)
    assert choose("FC(F)(F)CCCCCC")[0] == H.DLPNO


@pytest.mark.parametrize(
    "smiles, multiplicity, why",
    [
        ("C=CCCCCCCC", 1, "pi system"),
        ("c1ccc2ccccc2c1", 1, "pi system"),
        ("CCCCCCCC(=O)OC", 1, "pi system"),
        ("FC(F)(F)C(F)(F)CCCCC", 1, "fluorinated carbons"),
        ("[CH2]CCCCCCC", 2, "open-shell"),
        ("CCCCCCCCS", 1, "the number of cc-pVTZ-F12 basis functions of element 16"),
    ],
)
def test_no_high_result(smiles, multiplicity, why):
    level, reason = choose(smiles, multiplicity)
    assert level is None
    assert why in reason
