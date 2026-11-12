# -*- coding: utf-8 -*-

"""Matching a target to ATcT across resonance forms.

ATcT and the structure's SMILES may write one species as different resonance
forms: SiO is ``[Si]=O`` in ATcT but ``[O+]#[Si-]`` from a structure, and OBrO is
``[O]Br[O]`` against ``[O][Br+][O-]``. The canonical SMILES differ, so the
ordinary lookup misses. The key here ignores what differs between resonance
forms -- bond orders, formal charges and radical centres -- and keeps the atoms,
their connections (hydrogens included) and the total charge, so isomers such as
HCN and HNC still differ. A match is used only when exactly one ATcT species has
the key.
"""

import logging

logger = logging.getLogger(__name__)

#: {ATcT directory: {key: {flat SMILES, ...}}}, built on first use
_KEYS = {}


def resonance_key(smiles):
    """(the canonical SMILES of the bare graph, total charge), or None if the
    SMILES cannot be read."""
    from rdkit import Chem
    from seamm_thermochemistry.smiles import mol_from_smiles

    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    mol = Chem.RWMol(mol)
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol)
    graph = Chem.RWMol()
    for atom in mol.GetAtoms():
        bare = Chem.Atom(atom.GetAtomicNum())
        bare.SetNoImplicit(True)
        graph.AddAtom(bare)
    for bond in mol.GetBonds():
        graph.AddBond(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE
        )
    graph = graph.GetMol()
    graph.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(graph)
    return Chem.MolToSmiles(graph), charge


def _keys(tables):
    directory = str(tables.directory)
    if directory not in _KEYS:
        keys = {}
        for flat in tables.index["flat"]:
            try:
                key = resonance_key(flat)
            except Exception:  # an entry RDKit cannot handle this way
                key = None
            if key is not None:
                keys.setdefault(key, set()).add(flat)
        _KEYS[directory] = keys
    return _KEYS[directory]


def resonance_matches(tables, smiles, charge):
    """The ATcT SMILES (flat, canonical) of the gas-phase species with this
    charge that are the same species as `smiles` up to resonance."""
    key = resonance_key(smiles)
    if key is None:
        return []
    return sorted(
        flat
        for flat in _keys(tables).get(key, ())
        if tables.candidates(flat, charge=charge)
    )


def lookup_target(tables, smiles, charge, multiplicity):
    """The ATcT entry of the target, also when ATcT writes another resonance form.

    Returns
    -------
    (ATcTEntry or None, str or None)
        The entry, and the ATcT SMILES it was matched through when that is not the
        target's own (None for a direct match or no match).
    """
    from seamm_thermochemistry.atct import AmbiguousEntry, MissingEntry

    try:
        return tables.lookup(smiles, charge=charge, multiplicity=multiplicity), None
    except AmbiguousEntry:
        return None, None
    except MissingEntry:
        pass
    matches = resonance_matches(tables, smiles, charge)
    if len(matches) != 1:
        if len(matches) > 1:
            logger.info(
                f"{smiles} matches several ATcT species up to resonance: "
                + ", ".join(matches)
            )
        return None, None
    try:
        entry = tables.lookup(matches[0], charge=charge, multiplicity=multiplicity)
    except (MissingEntry, AmbiguousEntry):
        return None, None
    return entry, matches[0]
