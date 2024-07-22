"""
A collection of commonly used file utility functions
used when generating features from ColabFold complex
representations.
"""

import re
import os


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


def find_pdb_files(batch_path: str, stats_file: str) -> list[str]:
    """
    Crawl through a batch directory and return the PDB file
    that corresponds to the best model predicted by ColabFold.

    Parameters
    ----------
    batch_path : str
        The path to the batch directory you wish to extract
        PDB files from
    stats_file : str
        The path to the stats file that shows, for each complex,
        the number of the predicted model that is most confident

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
        pattern = rf"^{symbol}.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_\d+_seed_\d+\.pdb$"
        pdb_list.append(
            [file for file in os.listdir(batch_path) if re.match(pattern, file)][0]
        )
    return pdb_list
