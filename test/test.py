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
    expected_header = '>ENSG00000002822|ENST00000421113'
    expected_sequence = 'MRAILGSYDSELTPAEYSPQLTRRMREAEDMVQKVHSHSAE'
    actual_headers, actual_sequences = fasta_one_to_many.process_file(
        'test/test_data/single_gene.txt')
    assert (actual_headers[6] == expected_header and
            actual_sequences[6] == expected_sequence)
    