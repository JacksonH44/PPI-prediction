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


def apply_residue_mask(residue_sas: list[float], mask: list[bool]) -> tuple[int, int]:
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
    tuple[int, int]
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


def surface_area_stats(complex: str, complex_dir: str, monomer_dir: str):
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

    Returns
    -------
    list[str]
        A 12 element list consisting of avg, max, and min change in surface area
        for the interacting and non-interacting sites of both proteins in the
        complex
    """
    multimer_pdb_files = find_pdb_files(complex_dir)
    multimer_pdb = [file for file in multimer_pdb_files if complex in file]
    # Ensure that exactly 1 PDB file is found for the complex
    assert len(multimer_pdb) == 1, f'Unable to find 1 file for {complex} in {complex_dir}'
    multimer_pdb = multimer_pdb[0]
    logging.debug(f'Found multimer PDB file: {multimer_pdb}')
    # Find all potential monomer PDB files that might be needed
    monomers = os.listdir(monomer_dir)
    monomer_pdb_files = [pdb_file for m in monomers for pdb_file in find_pdb_files(os.path.join(monomer_dir, m))]

    symbol = multimer_pdb.split("/")[-1].split(".")[0]
    split = find_length_split(symbol, complex_dir)
    if split:
        file_path = os.path.join(complex_dir, multimer_pdb)
        seq_length_1, seq_length_2 = split[0], split[1]
        logging.debug(f"Generating contact map for {symbol}...")
        cmap_df = get_contact_map(file_path)
        logging.debug(f"Finding interaction site for {symbol}...")
        mask_map = find_interaction_site(
                symbol, cmap_df, seq_length_1, seq_length_2
        )
        logging.debug(f"Calculating multimer surface area structure for {symbol}...")
        multimer_struct = calculate_surface_areas(file_path)
        features = []
        # Perform the same commands for each of the two proteins in the complex
        for i in [0, 1]:
            sym = symbol.split("_")[i]
            # Find the corresponding folded monomer PDB file from ColabFold for the protein in the complex
            monomer_file = [pdb_file for pdb_file in monomer_pdb_files if sym in pdb_file]
            assert len(monomer_file) == 1, f'Could not find specific monomer file for {sym}.\nFound:\n{monomer_file}'
            monomer_symbol = monomer_file[0].split('.')[0]

            logging.debug(f'Calculating surface area structure for {monomer_symbol}...')
            monomer_struct = calculate_surface_areas(os.path.join(monomer_dir, monomer_symbol, monomer_file[0]))
            logging.debug(f"Extracting residue-level surface area for {sym}...")
            multimer_chain_sa = extract_residues(multimer_struct[0]["A" if i == 0 else "B"])
            monomer_chain_sa = extract_residues(monomer_struct[0]['A'])
            logging.debug(
                f"Splitting surface areas into interaction and non-interaction for {sym}..."
            )
            multimer_interaction, multimer_non_interaction = apply_residue_mask(
                multimer_chain_sa, mask_map[sym]
            )
            monomer_interaction, monomer_non_interaction = apply_residue_mask(
                monomer_chain_sa, mask_map[sym]
            )
            i_sa_avg, i_sa_max, i_sa_min = calculate_sa_metrics(monomer_interaction, multimer_interaction)
            ni_sa_avg, ni_sa_max, ni_sa_min = calculate_sa_metrics(monomer_non_interaction, multimer_non_interaction)
            features += [i_sa_avg, i_sa_max, i_sa_min, ni_sa_avg, ni_sa_max, ni_sa_min]
            
        driver = symbol.split('_')[0]
        interactor = symbol.split('_')[1]
        logging.debug(f'{symbol} stats:')
        logging.debug(f'{driver}_i_sa_avg,{driver}_i_sa_max,{driver}_i_sa_min,{driver}_ni_sa_avg,{driver}_ni_sa_max,{driver}_ni_sa_min,{interactor}_i_sa_avg,{interactor}_i_sa_max,{interactor}_i_sa_min,{interactor}_ni_sa_avg,{interactor}_ni_sa_max,{interactor}_ni_sa_min')
        features_str = [str(f) for f in features]
        logging.debug(','.join(features_str))
        return features


def find_delta_surface_areas(batch_dir: str, monomer_dir: str) -> None:
    """
    Find the change in surface area for the interaction
    and non-interaction sites for a batch of ColabFold outputs.

    Parameters
    ----------
    batch_dir : str
        The path to the file storing complexes of interest
    monomer_dir : str
        The path to the file storing monomer complexes of interest
    """
    multimer_pdb_files = find_pdb_files(batch_dir)
    monomers = os.listdir(monomer_dir)
    monomer_pdb_files = [pdb_file for m in monomers for pdb_file in find_pdb_files(os.path.join(monomer_dir, m))]
    for file in multimer_pdb_files:
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
            logging.debug(f"Calculating multimer surface area structure for {symbol}...")
            multimer_struct = calculate_surface_areas(file_path)
            for i in [0, 1]:
                sym = symbol.split("_")[i]
                monomer_file = [pdb_file for pdb_file in monomer_pdb_files if sym in pdb_file]
                assert len(monomer_file) == 1, f'Could not find specific monomer file for {sym}.\nFound:\n{monomer_file}'
                monomer_symbol = monomer_file[0].split('.')[0]
                logging.debug(f'Calculating surface area structure for {monomer_symbol}...')
                monomer_struct = calculate_surface_areas(os.path.join(monomer_dir, monomer_symbol, monomer_file[0]))
                logging.debug(f"Extracting residue-level surface area for {sym}...")
                multimer_chain_sa = extract_residues(multimer_struct[0]["A" if i == 0 else "B"])
                monomer_chain_sa = extract_residues(monomer_struct[0]['A'])
                logging.debug(
                    f"Splitting surface areas into interaction and non-interaction for {sym}..."
                )
                multimer_interaction, multimer_non_interaction = apply_residue_mask(
                    multimer_chain_sa, mask_map[sym]
                )
                monomer_interaction, monomer_non_interaction = apply_residue_mask(
                    monomer_chain_sa, mask_map[sym]
                )
                i_sa_avg, i_sa_max, i_sa_min = calculate_sa_metrics(monomer_interaction, multimer_interaction)
                ni_sa_avg, ni_sa_max, ni_sa_min = calculate_sa_metrics(monomer_non_interaction, multimer_non_interaction)
                logging.debug(f'Interaction site metrics - avg SA: {i_sa_avg}, max SA: {i_sa_max}, min SA: {i_sa_min}')
                logging.debug(f'Non-interaction site metrics - avg SA: {ni_sa_avg}, max SA: {ni_sa_max}, min SA: {ni_sa_min}')


def parse_command_line():  # pragma: no cover
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "batch_dir",
        type=str,
        help="The path to the batch directory holding files you want to get the surface area for",
    )
    parser.add_argument(
        '-m',
        '--monomer_dir',
        default=None,
        type=str,
        help='The path to the directory that holds ColabFold folded monomers'
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
    monomer_dir = args.monomer_dir if args.monomer_dir is not None else os.path.join(os.path.dirname(args.batch_dir), 'monomer')
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    # find_delta_surface_areas(args.batch_dir, monomer_dir)
    surface_area_stats('CDKN2A_CYCS', 'tests/test_data/colabfold/0', 'tests/test_data/colabfold/monomer')


if __name__ == "__main__":  # pragma: no cover
    main()
