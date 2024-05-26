import pytest
import os
import sys

bin_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../bin'))
sys.path.append(bin_dir)

import fasta_one_to_many


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
    