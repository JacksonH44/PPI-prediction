"""
A collection of functions to extract the change
in surface area from a complex.
"""

import argparse
import logging
import os
import sys
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
) -> tuple[float, float, float]:
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
    metrics : tuple[float, float, float]
        A triple of average change in SA, max change in SA, and min change in SA
    """
    # Check residues are the same length
    assert len(monomer_residues) == len(multimer_residues)
    deltas = []
    n = len(monomer_residues)
    for i in range(n):
        # Subtract monomer from multimer so positive delta
        # represents an increase in SA, and negative represents decrease
        deltas.append(multimer_residues[i] - monomer_residues[i])
    avg = round(sum(deltas) / n, 4)
    d_max = round(max(deltas), 4)
    d_min = round(min(deltas), 4)
    return (avg, d_max, d_min)


def split_residues_from_mask():
    """
    Get a list of residues and a mask as input, and return two lists of residues -
    one for the interaction site, and one for the non-interaction site
    """
    # TODO: finish function


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


def find_delta_surface_areas(batch_dir: str) -> None:
    """Find the change in surface area for the interaction
    and non-interaction sites for a batch of ColabFold outputs."""
    pdb_files = find_pdb_files(batch_dir)
    for file in pdb_files:
        symbol = file.split("/")[-1].split(".")[0]
        split = find_length_split(symbol, batch_dir)
        if split:
            file_path = os.path.join(batch_dir, file)
            seq_length_1, seq_length_2 = split[0], split[1]
            logging.debug(f"Generating contact map for {symbol}...")
            cmap_df = get_contact_map(file_path)
            logging.debug(f"Finding interaction site for {symbol}...")
            mask_map = find_interaction_site(
                symbol, cmap_df, seq_length_1, seq_length_2
            )
            logging.debug(f"Calculating surface area structure for {symbol}...")
            struct = calculate_surface_areas(file_path)
            logging.debug(
                f'Extracting residue-level surface area for {symbol.split("_")[0]}...'
            )
            chain_a_sa = extract_residues(struct[0]["A"])
            logging.debug(
                f'Extracting residue-level surface area for {symbol.split("_")[1]}...'
            )
            chain_b_sa = extract_residues(struct[0]["B"])


def parse_command_line():  # pragma: no cover
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "batch_dir",
        type=str,
        help="The path to the batch directory holding files you want to get the surface area for",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Change logging level from default level to noisiest level",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        type=str,
        default=None,
        help="Specify the name of the log file",
    )
    args = parser.parse_args()
    return args


def main():  # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    logfile = (
        args.logfile
        if args.logfile is not None
        else os.path.join(os.getcwd(), "logs/surface_area.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    find_delta_surface_areas(args.batch_dir)


if __name__ == "__main__":  # pragma: no cover
    main()
