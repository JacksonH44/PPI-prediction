"""
A test file for functions in the src/features directory.
"""

import pytest

import pandas as pd

from src.features.run_colabfold import create_observations, find_msa, prep_msas


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
