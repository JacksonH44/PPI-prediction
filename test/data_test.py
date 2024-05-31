import pytest
import os
import sys

import pandas as pd

bin_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/data'))
sys.path.append(bin_dir)

import fasta_one_to_many
import create_protein_triplets


def test_observation_creation1():
    """Test that the observation dataframe creation function works, includes the 
    correct labels for each triplet. Test 1/2."""
    expected_data = {
        'ref_ID': ['ACTN4_1', 'ACTN4_1', 'AKT1_1', 'AKT1_1', 'AKT1_1', 'BAG1_1'],
        'alt_ID': ['ACTN4_4', 'ACTN4_4', 'AKT1_2', 'AKT1_2', 'AKT1_2', 'BAG1_2'],
        'bait_ID': ['TRIM23', 'MYOZ2', 'TCL1A', 'TMCC2', 'MTUS2', 'HSPA8'],
        'perturbation': [True, True, False, True, True, True]
    }
    expected_df = pd.DataFrame(data=expected_data)
    actual_df = create_protein_triplets.find_triplets('test/test_data/test_ppis1.xlsx')
    assert expected_df.equals(actual_df)


def test_observation_creation2():
    """Test 2/2."""
    expected_data = {
        'ref_ID': ['BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1', 'BCL2L1_1',
                   'BCL2L1_1', 'BCL2L1_1', 'BTC_1', 'SERPING1_1', 'SERPING1_1'],
        'alt_ID': ['BCL2L1_2', 'BCL2L1_2', 'BCL2L1_2', 'BCL2L1_2', 'BCL2L1_3', 'BCL2L1_3',
                   'BCL2L1_3', 'BCL2L1_3', 'BTC_2', 'SERPING1_3', 'SERPING1_4'],
        'bait_ID': ['BAD', 'BIK', 'VAC14', 'BMF', 'BAD', 'BIK', 'VAC14', 'BMF', 'GMPPA',
                    'SLC30A2', 'SLC30A2'],
        'perturbation': [True, True, False, True, True, False, True, False, True, True, False]
    }
    expected_df = pd.DataFrame(data=expected_data)
    actual_df = create_protein_triplets.find_triplets('test/test_data/test_ppis2.xlsx')
    assert expected_df.equals(actual_df)


def test_filenotfounderror_handling():
    """Error handling test for a file that doesn't exist."""
    with pytest.raises(FileNotFoundError):
        actual_fname, actual_seq = fasta_one_to_many.process_file('../data/nonexistent.txt')


def test_file_separation():
    """Test that FASTA files are correctly separated 
    into the header and the protein sequence."""
    expected_header = '>ENSG00000002822|ENST00000455998'
    expected_sequence = '''MEDLGENTMVLSTLRSLNNFISQRVEGGSGLDISTSAPGSLQMQYQQSMQREVDRNQELL
TRIRQLQEREAGAEEKMQEQLERNRQCQQNLDAASKRLREKEDSLAQAGETINALKGRIS
ELQWSVMDQEMRVKRLESEKQELQ'''
    actual_headers, actual_sequences = fasta_one_to_many.process_file(
        'test/test_data/test_gene1.txt')
    assert (actual_headers[1] == expected_header and
            actual_sequences[1] == expected_sequence)
    

def test_file_creation():
    """Test the creation of individual FASTA files 
    from header sequence lists."""
    outfolder_name = os.path.join(os.getcwd(), 'test', 'test_data', 
                                  'test_gene1')
    fname = os.path.join(outfolder_name, 'ENSG00000002822_ENST00000421113.fasta')
    expected_content = '''>ENSG00000002822|ENST00000421113
MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE'''
    header = ['>ENSG00000002822|ENST00000421113']
    sequence = ['MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE']
    fasta_one_to_many._create_files(header, sequence, outfolder_name)
    with open(fname, 'r') as file:
        actual_content = file.read()
    assert actual_content == expected_content
    