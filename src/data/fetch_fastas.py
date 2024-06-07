"""
Fetch amino acid sequences in FASTA file format
from accession numbers.
"""


import argparse
import logging
import os

import pandas as pd
import requests

from core import config as cfg

API_KEY = cfg.NCBI_API_KEY
BASE_URL = cfg.NCBI_BASE_URL


def fetch_fastas(input_file, output_folder):
    """
    Make an NCBI API request for the amino acid sequence in
    a FASTA file for each EntrezGene accession number in the list.

    Parameters
    ----------
    reference file : str
        Absolute path to the input file that is an excel file with headers Isoform_ID
        and GenBank_Accession. An example of such file is 'ppis.xlsx'. 
    output_folder : str
        Absolute path to the output folder where all individual FASTA files
        will be written to
    """
    id_df = pd.read_excel(input_file, sheet_name='2A-Isoforms tested in Y2H', usecols=['Isoform_ID', 'GenBank_Accession'])
    accession_numbers = id_df['GenBank_Accession'].tolist()
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
            isoform_id = (id_df[id_df['GenBank_Accession'] == accession]['Isoform_ID'].values)[0]
            output_file = os.path.join(output_folder, f'{isoform_id}.fasta')
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
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/fetch_fastas.log")
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    fetch_fastas(args.infile, args.outfolder)


if __name__ == '__main__':
    main()
