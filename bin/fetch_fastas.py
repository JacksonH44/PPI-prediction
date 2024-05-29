"""
Fetch amino acid sequences in FASTA file format
from accession numbers.
"""


import argparse
import logging
import os

from dotenv import load_dotenv
import requests

load_dotenv()

API_KEY = os.getenv('API_KEY')
# Base URL for NCBI efetch
BASE_URL = os.getenv('EFETCH_BASE_URL')


def fetch_fastas(accession_numbers, output_folder):
    """Make an NCBI API request for the amino acid sequence in
    a FASTA file for each EntrezGene accession number in the list"""
    for accession in accession_numbers:
        params = {
            'db': 'nuccore',
            'id': accession,
            'rettype': 'fasta_cds_aa',
            'retmode': 'text',
            'api_key': API_KEY
        }
        response = requests.get(BASE_URL, params=params)
        logging.debug(f'Request for transcript from {accession} had API response {response.status_code}')
        if response.status_code == 200:
            fasta_sequence = response.text
            output_file = os.path.join(output_folder, f'{accession}.fasta')
            with open(output_file, "w") as file:
                file.write(fasta_sequence)
        else:
            logging.warning(f"Failed to retrieve sequence for {accession}")


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Input file name. FASTA files will be saved as their isoform names')
    parser.add_argument('outfolder', type=str,
                        help='Output folder path (requires absolute path)')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Change logging level from default level to noisiest level')
    parser.add_argument('-l', '--logfile',
                        type=str, default=None,
                        help='Specify the name of the log file')
    args = parser.parse_args()
    return args


def main():
    """Run the command line program."""
    args = parse_command_line()
    logfile = args.logfile if args.logfile is not None else f'{os.path.join(os.getcwd(), "logs/fetch_fastas.log")}'
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    fetch_fastas(['KU177872'], args.outfolder)


if __name__ == '__main__':
    main()
