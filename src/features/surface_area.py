"""
Compute the average, max, and min surface area across multiple
multimer models.
"""

import os
import re
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_pdb_files
from src.features.surface_area_calculator import SurfaceAreaCalculator


def surface_area_stats(
    complex: str, complex_dir: str, monomer_dir: str, num_models: int = 1, feature_type: str
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
    feature_type : str
        The type of feature you'd like to calculate

    Returns
    -------
    list[str]
        A 12 element list in the following order (i: interaction site, ni: non-interaction site,
        sa: surface area, 1: cancer driver gene, 2: interacting partner):
        [min_i_sa_1, min_i_sa_2, max_i_sa_1, max_i_sa_2, avg_i_sa_1, avg_i_sa_2,
        min_ni_sa_1, min_ni_sa_2, max_ni_sa_1, max_ni_sa_1, avg_ni_sa_1, avg_ni_sa_2]

    Usage
    -----
    surface_area_stats(
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
        absolute_multimer_pdb_path = os.path.join(complex_dir, multimer_pdb_file)
        for i in [0, 1]:
            symbol = complex.split("_")[i]
            monomer_pattern = rf".*_{symbol}.msa_unrelaxed_rank_001_alphafold2_ptm_model_\d+_seed_\d+\.pdb"
            monomer_res = [
                pdb_file
                for pdb_file in monomer_pdb_files
                if re.match(monomer_pattern, pdb_file.split("/")[-1])
            ]
            if len(monomer_res) == 1:
                monomer_file = monomer_res[0]
            else:
                monomer_res = [
                    pdb_file
                    for pdb_file in monomer_pdb_files
                    if symbol in pdb_file
                ][0]
            absolute_monomer_pdb_path = os.path.join(monomer_dir, monomer_file)
            sac = SurfaceAreaCalculator(absolute_multimer_pdb_path, absolute_monomer_pdb_path)
            sac.calculate_residue_metrics()
            interaction_delta, non_interaction_delta = sac.calculate_delta_metrics()
            # Calculate min, max, and avg delta SA for interaction and non-interaction site
            min_i_sa[i] = min(min_i_sa[i], interaction_delta)
            max_i_sa[i] = max(max_i_sa[i], interaction_delta)
            avg_i_sa[i] += interaction_delta
            min_ni_sa[i] = min(min_ni_sa[i], non_interaction_delta)
            max_ni_sa[i] = max(max_ni_sa[i], non_interaction_delta)
            avg_ni_sa[i] += non_interaction_delta

    avg_i_sa[0] = round(avg_i_sa[0] / num_models if avg_i_sa[0] != 0 else 0, 4)
    avg_i_sa[1] = round(avg_i_sa[1] / num_models if avg_i_sa[1] != 0 else 0, 4)
    avg_ni_sa[0] = round(avg_ni_sa[0] / num_models if avg_ni_sa[0] != 0 else 0, 4)
    avg_ni_sa[1] = round(avg_ni_sa[1] / num_models if avg_ni_sa[1] != 0 else 0, 4)
    features = []
    features += min_i_sa
    features += max_i_sa
    features += avg_i_sa
    features += min_ni_sa
    features += max_ni_sa
    features += avg_ni_sa
    return features


if __name__ == '__main__':
    features = surface_area_stats(
        'CDKN2A_CYCS',
        'tests/test_data/colabfold/0',
        'tests/test_data/colabfold/monomer',
        5
    )
    print(features)
