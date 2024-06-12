import pytest
import os

import aiohttp
import pandas as pd

from src.data.create_protein_triplets import find_triplets
from src.data.data_processing import remove_ground_truth_data
from src.data.fasta_one_to_many import _create_files, process_file
from src.data.generate_positive_dataset import chunk_input_genes, get_interactors, parse_input_genes


def test_remove_ground_truth_data_success():
    """Test the pruning of data in a dataset that 
    also appears in a ground truth dataset."""
    unpruned_ppis = ["NAT2_TEST", "TEST_ABCA3", "TEST-1_TEST-2", "SERPINA3_ACAA1"]
    expected_output = ["TEST-1_TEST-2"]
    actual_output = remove_ground_truth_data(
        unpruned_ppis, 
        'test/test_data/isoform_sequences_test.xlsx',
        '1A-Gene List',
        'Gene_Symbol',
        'data/interim/triplets.csv'
    )
    assert actual_output == expected_output


def test_chunk_input_genes_success():
    """Test the chunking of a list of genes into chunks of 
    genes (makes it easier for API calls)."""
    input_genes = ["FANCE", "BRCA1", "ARID1A", "BCL10", "ERBB4"]
    expected_chunked_genes = [
        ["FANCE", "BRCA1"],
        ["ARID1A", "BCL10"],
        ["ERBB4"]
    ]
    actual_chunked_genes = chunk_input_genes(input_genes, 2)
    assert actual_chunked_genes == expected_chunked_genes


@pytest.mark.parametrize("infile", [('test/test_data/cancer_driver_gene_list_test.csv')])
def test_parse_input_genes_success(infile):
    """Test the generation of reference gene list for PPI pairs."""
    expected_output = ["CHD4", "CHEK2", "CIC", "CIITA", "CLIP1", "CLTC", 
                       "CLTCL1", "CNBP", "CNOT3", "CNTRL", "COL1A1", 
                       "COL2A1", "CREB1", "CREB3L1", "CREB3L2", "CREBBP",
                       "CRLF2", "CRTC1"]
    actual_output = parse_input_genes(infile)
    assert expected_output == actual_output


@pytest.mark.asyncio
@pytest.mark.parametrize("gene_list, expected_interactions, threshold_level",
    [
        (
            ["FANCE"],
            {"FANCE_FANCC", "FANCE_FANCA", "FANCE_FANCD2", "FANCE_FANCF", 
             "FANCE_FANCM", "FANCE_FANCG", "FANCE_HES1", "FANCE_APITD1",
             "FANCE_STRA13"},
             2
        ),
        (
            ["ASXL1", "SS18"],
            {"ASXL1_BAP1", "ASXL1_FOXK1", "ASXL1_FOXK2", "ASXL1_HCFC1",
             "ASXL1_HIST1H1C", "ASXL1_AKT1", "SS18_SMARCA2", "SS18_SMARCB1",
             "SS18_RBM14"},
            3
        ),
        (
            ["TMPRSS2"],
            set(),
            2
        )
    ]
)
async def test_get_interactors_success(gene_list, expected_interactions, threshold_level):
    """Test creation of PPI interactions."""
    async with aiohttp.ClientSession() as session:
        ppis = await get_interactors(session, gene_list, threshold_level)
    actual_interactions = set(ppis)
    assert actual_interactions == expected_interactions


@pytest.mark.parametrize("test_file_path, expected_data, positive",
    [
        (
            'test/test_data/ppis_test1.xlsx', 
            {
                'ref_ID': ['ACTN4_1', 'ACTN4_1', 'AKT1_1', 'AKT1_1', 'AKT1_1', 'BAG1_1'],
                'alt_ID': ['ACTN4_4', 'ACTN4_4', 'AKT1_2', 'AKT1_2', 'AKT1_2', 'BAG1_2'],
                'bait_ID': ['TRIM23', 'MYOZ2', 'TCL1A', 'TMCC2', 'MTUS2', 'HSPA8'],
                'perturbation': [True, True, False, True, True, True]
            }, 
            False
        ),
        (
            'test/test_data/ppis_test2.xlsx', 
            {
                'ref_ID': ['BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1',
                     'BCL2L1_1', 'BTC_1'],
                'alt_ID': ['BCL2L1_2', 'BCL2L1_2', 'BCL2L1_2', 'BCL2L1_3', 'BCL2L1_3',
                     'BCL2L1_3', 'BTC_2'],
                'bait_ID': ['BAD', 'BIK', 'BMF', 'BAD', 'BIK', 'BMF', 'GMPPA'],
                'perturbation': [True, True, True, True, False, False, True]
            },     
            True
        ),
        (
            'test/test_data/ppis_test3.xlsx',
            {
                'ref_ID': ['CLCN2_1', 'CLCN2_1', 'CLCN2_1', 'CLCN2_1', 'CLCN2_1'],
                'alt_ID': ['CLCN2_2', 'CLCN2_3', 'CLCN2_4', 'CLCN2_4', 'CLCN2_5'],
                'bait_ID': ['UBQLN1_ORF2', 'FHL3', 'FHL3', 'UBQLN1_ORF2', 'FHL3'],
                'perturbation': [True, True, False, False, False]
            },
            False
        )
    ]
)
def test_find_triplets_success(test_file_path, expected_data, positive):
    """Test that the observation dataframe creation function works, includes the 
    correct labels for each triplet."""
    expected_df = pd.DataFrame(data=expected_data)
    actual_df = find_triplets(test_file_path, positive)
    assert expected_df.equals(actual_df)


def test_process_file_error():
    """Error handling test for a file that doesn't exist."""
    with pytest.raises(FileNotFoundError):
        actual_fname, actual_seq = process_file('../data/nonexistent.txt')


def test_process_file_success():
    """Test that FASTA files are correctly separated 
    into the header and the protein sequence."""
    expected_header = '>ENSG00000002822|ENST00000455998'
    expected_sequence = '''MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQREVDRNQELL
TRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRIS
ELQWSVMDQEMRVKRLESEKQELQ'''
    actual_headers, actual_sequences = process_file(
        'test/test_data/test_gene1.txt')
    assert (actual_headers[1] == expected_header and
            actual_sequences[1] == expected_sequence)
    

def test__create_files_success():
    """Test the creation of individual FASTA files 
    from header sequence lists."""
    outfolder_name = os.path.join(os.getcwd(), 'test', 'test_data', 
                                  'test_gene1')
    fname = os.path.join(outfolder_name, 'ENSG00000002822_ENST00000421113.fasta')
    expected_content = '''>ENSG00000002822|ENST00000421113
MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE'''
    header = ['>ENSG00000002822|ENST00000421113']
    sequence = ['MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE']
    _create_files(header, sequence, outfolder_name)
    with open(fname, 'r') as file:
        actual_content = file.read()
    assert actual_content == expected_content
    