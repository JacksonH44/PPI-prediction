"""
A collection of functions to calculate local
energetic frustration in complexes.
"""


import logging
import os
import subprocess
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_pdb_files
from src.features.interaction_site import (
    get_contact_map,
    find_interaction_site,
    find_length_split
)


def calculate_frstindex(pdb_file: str, chain: str) -> float:
    """
    Calculate the average frstIndex across all residues of a PDB file.

    Parameters
    ----------
    pdb_file : str
        The absolute path to the pdb_file you want to compute frustration index for
    chain : str
        The chain you want to compute frustration index for. If the PDB file represents
        a complex, this will be either 'A' or 'B', if the PDB file represents a
        monomer, this will be 'A'

    Returns
    -------
    avg_frst_index : float
        The average frustration index across all residues in the chain, rounded to 4
        significant figures
    """
    subprocess.run([
        "Rscript",
        'src/features/compute_frustration.R',
        pdb_file,
        '/tmp',
        chain
    ], stdout=open(os.devnull, 'wb'))
    # Get just the file name from the pdb file path
    pdb_file_name = pdb_file.split('/')[-1].rstrip('.pdb')
    output_file = os.path.join(
        '/tmp',
        pdb_file_name + '_' + chain + '.done',
        'FrustrationData',
        pdb_file_name + '_' + chain + '.pdb_singleresidue'
    )
    frst_df = pd.read_csv(
        output_file,
        sep=' ',
        usecols=['FrstIndex']
    )
    avg_frst_index = round(frst_df['FrstIndex'].sum() / frst_df.shape[0], 4)
    return avg_frst_index


def frustration_stats(complex: str, complex_dir: str, monomer_dir: str) -> list[str]:
    """
    Return the avg, max, and min change in FrstIndex for interaction and
    non-interaction site on both proteins in a complex.

    Parameters
    ----------
    complex : str
        The name of the complex to get the frustration stats for
    complex_dir : str
        The path to the bottom-most directory that contains the complex
    monomer_dir : str
        The path to the top-level monomer directory that contains all
        monomer ColabFold outputs

    Returns
    -------
    list[str]
        A 12 element list consisting of avg, max, and min change in FrstIndex
        for the interacting and non-interacting sites of both proteins in the
        complex

    Usage
    -----
    e.g., surface_area_stats('CDKN2A_CYCS', 'tests/test_data/colabfold/0', 'tests/test_data/colabfold/monomer')
    """
    multimer_pdb_files = find_pdb_files(complex_dir)
    multimer_pdb = [file for file in multimer_pdb_files if complex in file]
    # Ensure that exactly 1 PDB file is found for the complex
    assert (
        len(multimer_pdb) == 1
    ), f"Unable to find 1 file for {complex} in {complex_dir}"
    multimer_pdb_file = multimer_pdb[0]
    logging.debug(f"Found multimer PDB file: {multimer_pdb_file}")
    # Find all potential monomer PDB files that might be needed
    monomers = os.listdir(monomer_dir)
    monomer_pdb_files = [
        pdb_file
        for m in monomers
        for pdb_file in find_pdb_files(os.path.join(monomer_dir, m))
    ]

    symbol = multimer_pdb_file.split("/")[-1].split(".")[0]
    split = find_length_split(symbol, complex_dir)
    if split:
        file_path = os.path.join(complex_dir, multimer_pdb_file)
        seq_length_1, seq_length_2 = split[0], split[1]
        logging.debug(f"Generating contact map for {symbol}...")
        cmap_df = get_contact_map(file_path)
        logging.debug(f"Finding interaction site for {symbol}...")
        mask_map = find_interaction_site(symbol, cmap_df, seq_length_1, seq_length_2)
        logging.debug(f"Calculating multimer surface area structure for {symbol}...")


if __name__ == '__main__':
    avg_frst_index = calculate_frstindex(
        '/cluster/projects/kumargroup/jackson/colabfold/0/CDKN2A_CYCS.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000.pdb',
        'A'
    )
    print(avg_frst_index)
