"""
A test file for functions in the src/features directory.
"""

import pytest

import pandas as pd

from src.features.collect_colabfold_stats import get_colabfold_metrics
from src.features.file_utils import find_all_complexes, find_pdb_files
from src.features.run_colabfold import create_observations, find_msa, prep_msas
from src.features.find_stats import find_stats

# from src.features.interaction_site import (
#     apply_residue_mask,
#     find_interaction_site,
#     find_length_split,
# )


def test_surface_area_stats_success():
    """Test that the correct surface area stats are collected from
    a complex."""
    expected_result = [
        -49.3577,
        -59.8454,
        -30.1946,
        -35.837,
        -42.59048,
        -46.51482,
        -3.2862,
        -1.2802,
        -1.2838,
        -0.749,
        -2.34166,
        -1.0765,
    ]
    actual_result = find_stats(
        "CDKN2A_CYCS",
        "tests/test_data/colabfold/0",
        "tests/test_data/colabfold/monomer",
        "surface_area",
        5,
    )
    assert actual_result == expected_result


def test_get_colabfold_metrics_success():
    """Test that the correct pLDDT and ipTM scores are collected from
    a folded complex."""
    with open("tests/test_data/colabfold/0/log.txt", "r") as log:
        lines = log.readlines()
        lines = [line.split(" ", maxsplit=2)[2] for line in lines]
        lines = [line.rstrip("\n") for line in lines]
        expected_result = [
            "2",
            "79.5",
            "76.18",
            "79.5",
            "73.4",
            "0.501",
            "0.4894",
            "0.501",
            "0.476",
            "0.174",
            "0.1518",
            "0.174",
            "0.123",
        ]
        actual_result = get_colabfold_metrics("CDKN2A_TRAPPC2L", lines)
        assert expected_result == actual_result


# def test_apply_residue_mask_fail():
#     """Test that an AssertionError is raised when the length of the surface
#     areas is not equal to the length of the mask provided."""
#     with pytest.raises(AssertionError):
#         apply_residue_mask([4.7823, 1.1008, 0.8923, 2.1129], [True, False, False])


# def test_apply_residue_mask_success():
#     """Test that a tuple of interaction site surface areas and non-interaction site
#     surface areas is correctly returned."""
#     expected_result = ([4.7823, 2.1129], [1.1008, 0.8923])
#     actual_result = apply_residue_mask(
#         [4.7823, 1.1008, 0.8923, 2.1129], [True, False, False, True]
#     )
#     assert expected_result == actual_result


# def test_calculate_sa_metrics_fail():
#     """Test that the function correctly throws an error when two residue lists
#     are not the same length."""
#     with pytest.raises(AssertionError):
#         calculate_sa_metrics([3.23423, 4.2342], [8.23432])


# def test_calculate_sa_metrics_success():
#     """Test that the function correctly calculates the surface area metrics."""
#     expected_result = -0.9901
#     actual_result = calculate_sa_metrics(
#         [3.324, 8.3289, 2.234], [1.3425, 0.23423, 9.34]
#     )
#     assert expected_result == actual_result


# def test_find_interaction_site():
#     """Test that the correct residues are found as part of the interaction site."""
#     cmap_data = {
#         "residue_1": [1, 2, 2, 3, 5, 6, 7, 9],
#         "residue_2": [1, 2, 5, 7, 5, 3, 8, 8],
#     }
#     cmap_df = pd.DataFrame(data=cmap_data)
#     expected_result = {
#         "TEST1": [False, True, True, False],
#         "TEST2": [True, True, True, False, False],
#     }
#     actual_result = find_interaction_site("TEST1_TEST2", cmap_df, 4, 5)
#     assert actual_result == expected_result


# def test_find_length_split_fail():
#     """Test that when a symbol doesn't exist the error is handled elegantly."""
#     expected_result = None
#     actual_result = find_length_split("FANCE_FANCF", "tests/test_data/colabfold/0")
#     assert expected_result == actual_result


# def test_find_length_split_success():
#     """Verify that lengths are found accurately for complexes."""
#     expected_result = (164, 100)
#     actual_result = find_length_split("SRSF3_GGTA1", "tests/test_data/colabfold/0")
#     assert actual_result == expected_result


def test_find_pdb_files_success():
    """Test that find_pdb_files finds all the
    correct PDB files corresponding to the best models."""
    expected_result = set(
        [
            "CDKN2A_TRAPPC2L.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb",
            "CDKN2C_CD24.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb",
            "SRSF3_GGTA1.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb",
            "CDKN2A_CYCS.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_5_seed_000.pdb",
        ]
    )
    actual_result = set(find_pdb_files("tests/test_data/colabfold/0"))
    assert actual_result == expected_result


def test_find_all_complexes_success():
    """Test that finding complexes within a batch
    is correct."""
    expected_result = set(
        ["CDKN2A_TRAPPC2L", "CDKN2C_CD24", "SRSF3_GGTA1", "CDKN2A_CYCS"]
    )
    actual_result = set(find_all_complexes("tests/test_data/colabfold/0"))
    assert actual_result == expected_result


def test_create_observations_success():
    """Test that creating a dataframe with the
    required bounds works correctly."""
    data = {
        "symbol": ["TMA7", "NOP10", "SERP1", "FXYD2"],
        "length": [64, 64, 66, 66],
        "batch_number": [2, 2, 2, 2],
    }
    expected_result = pd.DataFrame(data=data)
    actual_result = create_observations("data/interim/sequence_lengths.csv", 2)
    actual_result = actual_result.reset_index(drop=True)
    expected_result = expected_result.reset_index(drop=True)
    assert expected_result.equals(actual_result)


def test_find_msa_success():
    """Test that a desired MSA is able to be found."""
    expected_result = "ENST00000586790_CCDC106.msa.a3m"
    actual_result = find_msa("CCDC106", "tests/test_data/msas")
    assert expected_result == actual_result


@pytest.mark.parametrize(
    "input, expected_output",
    [
        ("CCDC106", "tests/test_data/msas/ENST00000586790_CCDC106.msa.a3m"),
        ("ABI1_ENAH", "tests/test_data/msas/multimer/ABI1_ENAH.msa.a3m"),
    ],
)
def test_prep_msas_success(input, expected_output):
    """Test that MSA files are created and placed in the
    correct directory."""
    actual_output = prep_msas(input, "tests/test_data/msas")
    assert expected_output == actual_output
