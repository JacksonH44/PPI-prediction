"""
A test file for functions in the src/features directory.
"""

import pytest
import yaml

import pandas as pd

from src.features.collect_colabfold_stats import get_colabfold_metrics
from src.features.file_utils import find_all_complexes, find_pdb_files
from src.features.find_stats import find_stats
from src.features.run_colabfold import create_observations, find_msa, prep_msas
from src.features.surface_area_calculator import SurfaceAreaCalculator


multimer_pdb_path = (
    'tests/test_data/colabfold/0/CDKN2A_TRAPPC2L'
    '.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb'
)
monomer_pdb_path = (
    'tests/test_data/colabfold/monomer/ENST00000304494_CDKN2A'
    '.msa_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb'
)


def expected_result(var: str):
    """
    Return the expected result for a feature test.

    Parameters
    ----------
    var : str
        The name of the test case defined in feature_test.yml
        you wish to get the expected result for.

    Returns
    -------
    The expected result of the test
    """
    with open("tests/feature_test.yml", "r") as configFile:
        data = configFile.read()
    data = yaml.load(data, Loader=yaml.FullLoader)
    return data[var]["expected_output"]


def test_calculate_delta_metrics_success():
    """Test that the delta calculation (difference in monomer and multimer)
    is correct."""
    sac = SurfaceAreaCalculator(
        multimer_pdb_path,
        monomer_pdb_path
    )
    sac.calculate_residue_metrics()
    deltas = sac.calculate_delta_metrics()
    assert deltas == (-38.9749, -4.1072)


def test_surface_area_stats_success():
    """Test that the correct surface area stats are collected from
    a complex."""
    expected_result = [
        "-46.2103",
        "-56.1303",
        "-30.961",
        "-35.837",
        "-40.6794",
        "-44.2668",
        "-2.8622",
        "-1.9769",
        "-2.1979",
        "-0.3732",
        "-2.4745",
        "-1.2583",
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
    expected_result = [
        "31.65625",
        "76.17",
        "79.5144",
        "73.3834",
        "0.49",
        "0.5",
        "0.48",
        "0.15",
        "0.17",
        "0.12",
    ]
    actual_result = get_colabfold_metrics(
        "CDKN2A_TRAPPC2L", "tests/test_data/colabfold/0/"
    )
    assert expected_result == actual_result


def test_apply_residue_mask_success():
    """Test that a tuple of interaction site surface areas and non-interaction site
    surface areas is correctly returned."""
    sac = SurfaceAreaCalculator(
        multimer_pdb_path,
        monomer_pdb_path
    )
    sac.calculate_residue_metrics()
    expected_result = ([], [4.7823, 1.1008, 0.8923, 2.1129])
    actual_result = sac._apply_residue_mask([4.7823, 1.1008, 0.8923, 2.1129])
    assert expected_result == actual_result


@pytest.mark.parametrize(
    "result_type, expected_result_variable",
    [
        ("monomer", "test_calculate_residue_metrics_monomer"),
        ("multimer", "test_calculate_residue_metrics_multimer"),
    ],
)
def test_calculate_residue_metrics_success(result_type, expected_result_variable):
    """Test that monomer and multimer residue metrics are calculated correctly."""
    sac = SurfaceAreaCalculator(
        multimer_pdb_path,
        monomer_pdb_path
    )
    sac.calculate_residue_metrics()
    actual_output = sac._monomer_residue_metrics
    if result_type == "multimer":
        actual_output = sac._multimer_residue_metrics
    assert actual_output == expected_result(expected_result_variable)


def test_create_interaction_site_success():
    """Test that a FeatureCalculator object accurately creates a mask map."""
    sac = SurfaceAreaCalculator(
        multimer_pdb_path,
        monomer_pdb_path
    )
    assert sac._mask_map == expected_result("test_create_interaction_site")


def test_create_length_split_success():
    """Test that the proper length split is created upon the creation
    of a FeatureCalculator object."""
    sac = SurfaceAreaCalculator(
        multimer_pdb_path,
        monomer_pdb_path
    )
    assert sac._seq_1_length == 156 and sac._seq_2_length == 139


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
