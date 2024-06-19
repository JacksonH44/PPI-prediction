"""A file that gets all MANE transcripts Gencode IDs for 
a dataset."""

import argparse
import csv
import logging
import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from core import config as cfg
from src.data.data_processing import find_unique_genes


def find_mane_transcripts(unique_genes) -> tuple[pd.DataFrame, set[str]]:
    """Create a dataframe of gene names and Gencode transcript IDs"""
    mane_df = pd.read_csv(cfg.MANE_FILE)
    unique_genes_df = mane_df[mane_df['gene'].isin(unique_genes)]
    non_transcript_genes = unique_genes - set(unique_genes_df['gene'].to_list())
    return (unique_genes_df, non_transcript_genes)


def parse_command_line():  # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'files',
        nargs='+',
        type=str,
        help='A list of CSV files that you want Gencode IDs for'
    )
    parser.add_argument(
        '-t',
        '--transcript_file',
        type=str,
        default='data/interim/dataset_gencode_transcripts.csv',
        help='A list of CSV files that you want Gencode IDs for'
    )
    parser.add_argument(
        '-n',
        '--no_transcript_file',
        type=str,
        default='data/interim/dataset_no_gencode_transcripts.csv',
        help='A list of CSV files that you want Gencode IDs for'
    )
    parser.add_argument(
        "-l",
        "--logfile",
        type=str,
        default=None,
        help="Specify the name of the log file",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Change logging level from default level to noisiest level",
    )
    args = parser.parse_args()
    return args


def main():  # pragma: no cover
    """Run the main command line program."""
    args = parse_command_line()
    logfile = (
        args.logfile
        if args.logfile is not None
        else os.path.join(os.getcwd(), "logs/find_gencode_ids.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    unique_genes = find_unique_genes(args.files)
    transcripts_df, no_transcripts = find_mane_transcripts(unique_genes)
    logging.debug(f'Found {len(transcripts_df["gene"].to_list())} genes with a MANE transcript...')
    logging.debug(f'Found {len(no_transcripts)} genes without a MANE transcript...')
    logging.debug(f'Writing genes with a transcript to {args.transcript_file}...')
    transcripts_df.to_csv(args.transcript_file, header=True, index=False)
    logging.debug(f'Writing genes without a transcript to {args.no_transcript_file}...')
    with open(args.no_transcript_file, 'w') as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(['gene'])
        [writer.writerow([row]) for row in list(no_transcripts)]


if __name__ == '__main__':  # pragma: no cover
    main()