"""
A collection of functions that find the interaction site
of a protein complex.
"""

import logging
from typing import Optional

import pandas as pd
import numpy as np
from contact_map import ContactFrequency  # type: ignore
import mdtraj as md  # type: ignore


def apply_residue_mask(
    residues: list[float], mask: list[bool]
) -> tuple[list[float], list[float]]:
    """
    Get a list of residues and a mask as input, and return two lists of residues -
    one for the interaction site, and one for the non-interaction site

    Parameters
    ----------
    residues : list[float]
        List of all residue metrics for a sequence
    mask : list[bool]
        Interaction site mask for the sequence

    Returns
    -------
    tuple[list[float], list[float]]
        A tuple of residue-level metrics in the order (interaction_site, non_interaction_site)
    """
    assert len(residues) == len(mask)
    interaction_site = [residues[i] for i in range(len(residues)) if mask[i]]
    non_interaction_site = [residues[i] for i in range(len(residues)) if not mask[i]]
    return (interaction_site, non_interaction_site)


def find_interaction_site(
    complex_symbol: str,
    contact_map_df: pd.DataFrame,
    seq_length_1: int,
    seq_length_2: int,
) -> dict[str, list[bool]]:
    """
    Find the interaction sites for the two proteins present in the complex. It considers
    cases where one residue in the contact map is on one protein, while the other residue
    is on the other protein. Indices returned in the mask are relative to the single protein
    (e.g. if an interacting residue is residue #214 in the complex but #89 in the single
    protein, the number 89 will be returned).

    Parameters
    ----------
    complex_symbol : str
        The name of the complex to find interaction sites for
    contact_map_df : pd.DataFrame
        A pandas dataframe with column 1 being the first residue index of interaction,
        column 2 being the second residue index of interaction
    seq_length_1 : int
        The length of the first protein
    seq_length_2 : int
        The length of the second protein

    Returns
    -------
    mask_map : dict[str, list[bool]]
        A two-element hashmap that contains symbol name, interaction site mask pairs
        where the interaction site mask is a list of booleans with the same length as
        the symbol sequence. Values are true for residues that are part of the
        interaction site, and false for those that aren't.
    """
    protein_1, protein_2 = complex_symbol.split("_")[0], complex_symbol.split("_")[1]
    mask_1 = [False] * seq_length_1
    mask_2 = [False] * seq_length_2
    for row in contact_map_df.itertuples(index=False):
        residue_1, residue_2 = row[0], row[1]
        # NOTE: Indices of the masks are one less than the residues because
        # residues are 1-indexed and masks are 0-indexed
        # residue 1 is on protein 1 and residue 2 on protein 2
        if residue_1 <= seq_length_1 and residue_2 > seq_length_1:
            mask_1[residue_1 - 1] = True
            mask_2[residue_2 - seq_length_1 - 1] = True
        # residue 1 is on protein 2 and residue 2 on protein 1
        if residue_1 > seq_length_1 and residue_2 <= seq_length_1:
            mask_2[residue_1 - seq_length_1 - 1] = True
            mask_1[residue_2 - 1] = True
    mask_map = {}
    mask_map[protein_1] = mask_1
    mask_map[protein_2] = mask_2
    return mask_map


def find_length_split(symbol_name: str, batch_dir: str) -> Optional[tuple[int, int]]:
    """
    Find the sequence lengths for the two proteins in a protein
    complex.

    Parameters
    ----------
    symbol_name : str
        The name of complex to find lengths for
    batch_dir : str
        The path to the directory that the MSA is stored in

    Returns
    -------
    tuple[int, int], optional
        A tuple of the cancer driver sequence length, and the interactor
        sequence length
    """
    msa_file = f"{batch_dir}/{symbol_name}.msa.a3m"
    try:
        with open(msa_file, "r") as msa:
            firstLine = msa.readline().strip("\n")
            lengths = firstLine.split("\t")[0].lstrip("#").split(",")
            seq_1, seq_2 = int(lengths[0]), int(lengths[1])
            return (seq_1, seq_2)
    except FileNotFoundError:
        logging.error(f"Cannot find processed MSA for {symbol_name} in {batch_dir}")
        return None


def get_contact_map(pdb_file: str) -> pd.DataFrame:
    """
    Get the contact map for a PDB file. Used script
    from NourHanafi on Github.

    Parameters
    ----------
    pdb_file : str
        The path to the PDB file for a complex to obtain the
        contact map for

    Returns
    -------
    pd.DataFrame
        A dataframe containing the residue numbers that are
        in contact
    """
    traj = md.load(pdb_file)
    frame_contacts = ContactFrequency(traj)
    df = frame_contacts.residue_contacts.df

    df = df.fillna(0)
    df = df.where(np.triu(np.ones(df.shape)).astype(bool))
    df = df.stack().reset_index()
    df.columns = ["residue1", "residue2", "weight"]

    df["residue1"] += 1
    df["residue2"] += 1

    df = df.loc[df.residue1.ne(df.residue2)]  # keep non-self rows
    df = df[df["weight"] == 1.0]
    return df[["residue1", "residue2"]]
