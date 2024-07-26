"""
A collection of functions to extract the change
in surface area from a complex.
"""

import logging
import os
import sys
import re
from typing import Optional

from Bio.PDB import PDBParser, SASA  # type: ignore
from contact_map import ContactFrequency  # type: ignore
import mdtraj as md  # type: ignore
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_pdb_files


def calculate_sa_metrics(
    monomer_residues: list[float], multimer_residues: list[float]
) -> float:
    """
    Return the average, minimum, and maximum change in surface area (SA) between
    two sets of residues

    Parameters
    ----------
    monomer_residues : list[float]
        The SAs for the residues when they are in monomer form
    multimer_residues : list[float]
        The SAs for the residues when they are in multimer form

    Returns
    -------
    metrics : float
        A triple of average change in SA, max change in SA, and min change in SA
    """
    assert len(monomer_residues) == len(
        multimer_residues
    ), "The lengths of monomer and multimer residues should be the same"
    return round(
        (sum(multimer_residues) - sum(monomer_residues)) / len(monomer_residues), 4
    )


def apply_residue_mask(
    residue_sas: list[float], mask: list[bool]
) -> tuple[list[float], list[float]]:
    """
    Get a list of residues and a mask as input, and return two lists of residues -
    one for the interaction site, and one for the non-interaction site

    Parameters
    ----------
    residue_sas : list[float]
        List of all residue surface areas for a sequence
    mask : list[bool]
        Interaction site mask for the sequence

    Returns
    -------
    tuple[list[float], list[float]]
        A tuple of residue-level surface areas in the order (interaction_site, non_interaction_site)
    """
    assert len(residue_sas) == len(mask)
    interaction_site = [residue_sas[i] for i in range(len(residue_sas)) if mask[i]]
    non_interaction_site = [
        residue_sas[i] for i in range(len(residue_sas)) if not mask[i]
    ]
    return (interaction_site, non_interaction_site)


def extract_residues(chain) -> list[float]:
    """
    Return a list with the same length as the input sequence of surface areas,
    one entry for each sequence. The model index and chain should already be specified
    (e.g. usage: extract_residues(struct[0]['A']))

    Parameters
    ----------
    chain : Chain
        The input chain for which to extract residue-level surface area calculations from

    Returns
    -------
    surface_areas : list[float]
        A list of surface areas the same length as the number of residues in the chain rounded
        to 4 significant digits
    """
    surface_areas = [round(chain[i].sasa, 4) for i in range(1, len(chain) + 1)]
    return surface_areas


def calculate_surface_areas(pdb_path: str):
    """
    Calculates surface areas for each residue of a structure in a PDB file

    Parameters
    ----------
    pdb_path : str
        The path to the PDB file you want to calculate surface areas of

    Returns
    -------
    sr : Structure
        A structure object representing the residue-level computed surface areas
        for the PDB file
    """
    p = PDBParser(QUIET=1)
    symbol = pdb_path.split("/")[-1].split(".")[0]
    struct = p.get_structure(symbol, pdb_path)
    sr = SASA.ShrakeRupley()
    sr.compute(struct, level="R")
    return struct


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
    """Find the sequence lengths for the two proteins in a protein
    complex."""
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
    """Get the contact map for a PDB file. Used script
    from NourHanafi on Github."""
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


def surface_area_stats(
    complex: str, complex_dir: str, monomer_dir: str, num_models: int = 1
) -> list[float]:
    """
    Return the avg, max, and min change in surface area for interaction
    and non-interaction site on a complex.

    Parameters
    ----------
    complex : str
        The name of the complex to get stats for
    complex_dir : str
        The path to the bottom-most directory that contains the complex
    monomer_dir : str
        The path to the top-level monomer directory that contains all monomer
        Colabfold outputs
    num_models : int
        The number of AlphaFold models you'd like to take into account (defaults to 1,
        use 5 if you'd like to include all models)

    Returns
    -------
    list[str]
        A 12 element list in the following order (i: interaction site, ni: non-interaction site,
        sa: surface area, 1: cancer driver gene, 2: interacting partner):
        [min_i_sa_1, min_i_sa_2, max_i_sa_1, max_i_sa_2, avg_i_sa_1, avg_i_sa_2,
        min_ni_sa_1, min_ni_sa_2, max_ni_sa_1, max_ni_sa_1, avg_ni_sa_1, avg_ni_sa_2]

    Usage
    -----
    e.g., surface_area_stats(
        'CDKN2A_CYCS',
        'tests/test_data/colabfold/0',
        'tests/test_data/colabfold/monomer',
        5
    )
    """
    all_pdb_files = find_pdb_files(complex_dir, num_models)
    multimer_pdb = [file for file in all_pdb_files if complex in file]
    monomers = os.listdir(monomer_dir)
    monomer_pdb_files = [
        os.path.join(m, pdb_file)
        for m in monomers
        for pdb_file in find_pdb_files(os.path.join(monomer_dir, m), num_models)
    ]
    split = find_length_split(complex, complex_dir)
    if split is None:
        return [float("nan")] * 12

    min_i_sa = [float("inf"), float("inf")]
    max_i_sa = [float("-inf"), float("-inf")]
    avg_i_sa = [0.0, 0.0]
    min_ni_sa = [float("inf"), float("inf")]
    max_ni_sa = [float("-inf"), float("-inf")]
    avg_ni_sa = [0.0, 0.0]
    for rank in range(1, num_models + 1):
        multimer_pattern = rf"^{complex}.msa_unrelaxed_rank_00{rank}_alphafold2_multimer_v3_model_\d+_seed_\d+\.pdb$"
        multimer_pdb_file = [
            file for file in multimer_pdb if re.match(multimer_pattern, file)
        ][0]
        # generate contact map for this model
        file_path = os.path.join(complex_dir, multimer_pdb_file)
        cmap_df = get_contact_map(file_path)
        # find the interaction site for this complex based on the contact map
        mask_map = find_interaction_site(complex, cmap_df, split[0], split[1])
        # calculate the surface areas for the complex
        multimer_struct = calculate_surface_areas(file_path)
        # analyze monomers
        for i in [0, 1]:
            symbol = complex.split("_")[i]
            monomer_pattern = rf".*_{symbol}.msa_unrelaxed_rank_00{rank}_alphafold2_ptm_model_\d+_seed_\d+\.pdb"
            monomer_file = [
                pdb_file
                for pdb_file in monomer_pdb_files
                if re.match(monomer_pattern, pdb_file.split("/")[-1])
            ][0]
            # calculate surface areas for the monomer
            monomer_struct = calculate_surface_areas(
                os.path.join(monomer_dir, monomer_file)
            )
            # extract residues from multimer chain
            multimer_chain_sa = extract_residues(
                multimer_struct[0]["A" if i == 0 else "B"]
            )
            # extract residues from monomer chain
            monomer_chain_sa = extract_residues(monomer_struct[0]["A"])
            # divide multimer and monomer SAs into interaction and non-interaction
            multimer_interaction, multimer_non_interaction = apply_residue_mask(
                multimer_chain_sa, mask_map[symbol]
            )
            monomer_interaction, monomer_non_interaction = apply_residue_mask(
                monomer_chain_sa, mask_map[symbol]
            )
            interaction_delta = calculate_sa_metrics(
                monomer_interaction, multimer_interaction
            )
            non_interaction_delta = calculate_sa_metrics(
                monomer_non_interaction, multimer_non_interaction
            )

            # Calculate min, max, and avg delta SA for interaction and non-interaction site
            min_i_sa[i] = min(min_i_sa[i], interaction_delta)
            max_i_sa[i] = max(max_i_sa[i], interaction_delta)
            avg_i_sa[i] += interaction_delta
            min_ni_sa[i] = min(min_ni_sa[i], non_interaction_delta)
            max_ni_sa[i] = max(max_ni_sa[i], non_interaction_delta)
            avg_ni_sa[i] += non_interaction_delta

    avg_i_sa[0] = avg_i_sa[0] / num_models if avg_i_sa[0] != 0 else 0
    avg_i_sa[1] = avg_i_sa[1] / num_models if avg_i_sa[1] != 0 else 0
    avg_ni_sa[0] = avg_ni_sa[0] / num_models if avg_ni_sa[0] != 0 else 0
    avg_ni_sa[1] = avg_ni_sa[1] / num_models if avg_ni_sa[1] != 0 else 0
    features = []
    features += min_i_sa
    features += max_i_sa
    features += avg_i_sa
    features += min_ni_sa
    features += max_ni_sa
    features += avg_ni_sa
    return features
