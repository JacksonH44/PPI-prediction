"""
A test file for functions in the src/features directory.
"""

import pytest

import pandas as pd

from src.features.run_colabfold import create_observations, find_msa, prep_msas
from src.features.file_utils import find_all_complexes, find_pdb_files
from src.features.surface_area import find_length_split


def test_find_length_split_fail():
    """Test that when a symbol doesn't exist the error is handled elegantly."""
    expected_result = None
    actual_result = find_length_split('FANCE_FANCF', 'tests/test_data/colabfold/0')
    assert expected_result == actual_result


def test_find_length_split_success():
    """Verify that lengths are found accurately for complexes."""
    expected_result = (164, 100)
    actual_result = find_length_split('SRSF3_GGTA1', 'tests/test_data/colabfold/0')
    assert actual_result == expected_result


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
    print(actual_result)
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
