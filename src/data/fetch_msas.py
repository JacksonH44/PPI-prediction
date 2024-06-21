"""
Find all MSAs from UniClust30 DB for each protein
in the dataset.
"""

import argparse
from collections import defaultdict
import logging
import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from core import config as cfg
from src.data.data_processing import find_unique_genes


def find_canonical_transcript(genes) -> defaultdict[str, str]:
    """Find the canonical transcripts for the whole set of genes."""
    transcripts = defaultdict(str)
    mane_df = pd.read_csv(cfg.MANE_FILE, sep="\t", usecols=["symbol", "Ensembl_nuc"])
    for gene in genes:
        transcripts[gene] = mane_df[mane_df["symbol"] == gene]["Ensembl_nuc"].values[0]
        if transcripts[gene] == "":  # Symbol not in MANE summary file
            logging.warning(f"Did not find canonical transcript for {gene}")
    return transcripts


def parse_command_line():  # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        nargs="+",
        type=str,
        help="The input file(s) specifying proteins that you want to find MSAs for",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="data/processed/negative_ppis.csv",
        help="Path to output file to store negative PPI datset",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Change logging level from default level to noisiest level",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        type=str,
        default=None,
        help="Specify the name of the log file",
    )
    args = parser.parse_args()
    return args


def main():  # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    logfile = (
        args.logfile
        if args.logfile is not None
        else os.path.join(os.getcwd(), "logs/fetch_msas.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    all_genes = find_unique_genes(args.input_file)
    logging.debug(f"Found {len(all_genes)} unique genes...")
    transcripts = find_canonical_transcript(all_genes)
    # Query BioMart API for UniRef30 protein IDs

    # Use BioPython to find longest protein sequence of the duplicates

    # Fetch MSAs and store them somewhere


if __name__ == "__main__":  # pragma: no cover
    main()
