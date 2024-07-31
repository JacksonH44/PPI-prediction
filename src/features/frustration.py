"""
A collection of functions to calculate local
energetic frustration in complexes.
"""


import logging
import os
import re
import subprocess
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_pdb_files
from src.features.interaction_site import (
    apply_residue_mask,
    get_contact_map,
    find_interaction_site,
    find_length_split
)


def calculate_avg_frst_index_delta(multimer_residues: list[float], monomer_residues: list[float]) -> float:
    """
    Calculate the change in frustration index between a multimer
    and monomer complex.

    Parameters
    ----------
    multimer_residues : list[float]
        The list of frustration indices for the multimer residues
    monomer_residues : list[float]
        The list of frustration indices for the monomer residues

    Returns
    -------
    float
        The difference in average frustration index between the multimer
        and monomer
    """
    assert len(monomer_residues) == len(
        multimer_residues
    ), "The lengths of monomer and multimer residues should be the same"
    return round((sum(multimer_residues) - sum(monomer_residues)) / len(monomer_residues), 4)


def get_frstindex(pdb_file: str, chain: str) -> list[float]:
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
    list[float]
        A list of residue-level frustration index values
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
    return frst_df['FrstIndex'].to_list()


def calculate_deltas(
        symbol: str,
        multimer_file: str, 
        monomer_file: str, 
        partner_num: int,
        mask_map: dict[str, list[bool]]
) -> tuple[float, float]:
    """
    Calculate the interaction and non-interaction site deltas for a gene in a complex.

    Parameters
    ----------
    symbol : str
        Gene symbol
    multimer_file : str
        The full path to the multimer file
    monomer_file : str
        The full path to the monomer file
    partner_num : int
        The number partner the gene is (either 0 or 1)
    mask_map : dict[str, list[bool]]
        A map of gene symbols to its interaction site mask

    Returns
    -------
    tuple[float, float]
        The delta values for the interaction and non-interaction site
    """
    multimer_frst_index = get_frstindex(multimer_file, 'A' if partner_num == 0 else 'B')
    monomer_frst_index = get_frstindex(monomer_file, 'A')
    # divide multimer and monomer SAs into interaction and non-interaction
    multimer_interaction, multimer_non_interaction = apply_residue_mask(
        multimer_frst_index, mask_map[symbol]
    )
    monomer_interaction, monomer_non_interaction = apply_residue_mask(
        monomer_frst_index, mask_map[symbol]
    )
    interaction_delta = calculate_avg_frst_index_delta(
        monomer_interaction, multimer_interaction
    )
    non_interaction_delta = calculate_avg_frst_index_delta(
        monomer_non_interaction, multimer_non_interaction
    )
    return (interaction_delta, non_interaction_delta)


def frustration_stats(
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
    frustration_stats(
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
    
    min_i_f = [float("inf"), float("inf")]
    max_i_f = [float("-inf"), float("-inf")]
    avg_i_f = [0.0, 0.0]
    min_ni_f = [float("inf"), float("inf")]
    max_ni_f = [float("-inf"), float("-inf")]
    avg_ni_f = [0.0, 0.0]



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

        # analyze monomers
        for i in [0, 1]:
            symbol = complex.split("_")[i]
            # calculate multimer frustration index
            multimer_frst_index = get_frstindex(
                os.path.join(complex_dir, multimer_pdb_file),
                'A' if i == 0 else 'B'
            )
            monomer_pattern = rf".*_{symbol}.msa_unrelaxed_rank_001_alphafold2_ptm_model_\d+_seed_\d+\.pdb"
            monomer_file = [
                pdb_file
                for pdb_file in monomer_pdb_files
                if re.match(monomer_pattern, pdb_file.split("/")[-1])
            ][0]
            # calculate frustration index for the monomer
            monomer_frst_index = get_frstindex(os.path.join(monomer_dir, monomer_file), 'A')


if __name__ == '__main__':
    frustration_stats(
        'CDKN2A_CYCS',
        'tests/test_data/colabfold/0',
        'tests/test_data/colabfold/monomer',
        5
    )
