"""
A collection of commonly used file utility functions
used when generating features from ColabFold complex
representations.
"""

import re
import os

import pandas as pd


def find_all_complexes(batch_path: str) -> list[str]:
    """
    Find all the complexes that are in a batch outputted
    by ColabFold batch.

    Parameters
    ----------
    batch_path : str
        The path to the batch directory holding the complex
        predictions you wish to extract

    Returns
    -------
    complexes : list[str]
        A list of all the unique complex names in the batch
    """
    omit = set(["cite.bibtex", "config.json", "log.txt"])
    complex_set = set()
    for file in os.listdir(batch_path):
        name = file.split(".")[0]
        if file not in omit:
            complex_set.add(name)

    complexes = list(complex_set)
    return complexes


def find_pdb_files(batch_path: str, num_models: int = 1) -> list[str]:
    """
    Crawl through a batch directory and return the PDB file
    that corresponds to the best model predicted by ColabFold.

    Parameters
    ----------
    batch_path : str
        The path to the batch directory you wish to extract
        PDB files from
    num_models : int
        The number of models to get PDB files for, will get the best n
        models you specify

    Returns
    -------
    pdb_list : list[str]
        A list of file paths the size of the batch size for that batch
        of PDB files corresponding to the best model for each complex
        predicted by ColabFold
    """
    pdb_list = []
    complex_symbols = find_all_complexes(batch_path)
    for symbol in complex_symbols:
        base_pdb_pattern = rf"^{symbol}\.pdb$"
        base_pdbs = [
            file for file in os.listdir(batch_path) if re.match(base_pdb_pattern, file)
        ]
        if len(base_pdbs) != 0:
            pdb_list.append(base_pdbs[0])
        for model_no in range(1, num_models + 1):
            multimer_pattern = (
                rf"^{symbol}.msa_unrelaxed_rank_00{model_no}_alphafold2_multimer_v3_model_\d+_seed_\d+\.pdb$"
            )
            monomer_pattern = (
                rf"^{symbol}.msa_unrelaxed_rank_00{model_no}_alphafold2_ptm_model_\d+_seed_\d+\.pdb"
            )
            files_found = [
                file
                for file in os.listdir(batch_path)
                if re.match(multimer_pattern, file) or re.match(monomer_pattern, file)
            ]
            if len(files_found) != 0:
                pdb_list.append(files_found[0])
    return pdb_list


def combine_csv(data_dir: str, upper_bound: str) -> None:
    """
    Combine multiple CSVs into one CSV, assuming they all have
    the same headers. It also assumes that the CSVs are labeled
    colabfold_stats_x.csv where x is the batch of ColabFold outputs
    it represents.

    Parameters
    ----------
    data_dir : str
        The path to the directory holding all the intermediate CSVs
    upper_bound : str
        The largest number x for which the file colabfold_stats_x.csv
        exists for
    """
    dfs = []
    for i in range(upper_bound + 1):
        file = os.path.join(data_dir, f"colabfold_stats_{i}.csv")
        file_df = pd.read_csv(file)
        dfs.append(file_df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_csv("data/processed/colabfold_stats.csv", index=False)
