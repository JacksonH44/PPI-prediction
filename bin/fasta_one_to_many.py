"""
Filter and convert one big FASTA file into 
many smaller files for each transcript.
"""


import argparse
import os


def fasta_one_to_many(file_path):
    pass


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Input file name')
    parser.add_argument('-o', '--outfolder',
                        type=str, default=None,
                        help='Output folder path')
    args = parser.parse_args()
    return args


def main():
    """Run the command line program."""
    args = parse_command_line()
    outfolder = os.getcwd()

if __name__ == '__main__':
    main()