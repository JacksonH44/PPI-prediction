"""
Run ColabFold for protein pairs in a dataset.
"""

import argparse
import logging
import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.data.combine_msa import (
    extract_header_sequence_pairs,
    write_combined_a3m
)


def find_msa(gene_symbol, msa_dir) -> str:
    """Find the MSA for a given gene symbol in the MSA directory."""
    all_msas = os.listdir(msa_dir)
    msa_file = [msa for msa in all_msas if gene_symbol in msa]
    # There should only be one entry in msa_file
    if len(msa_file) == 1:
        return msa_file[0]
    else:
        logging.warning(f'Could not find MSA file for {gene_symbol}')
        return 'NA'


def prep_msas(symbol: str, msa_dir: str) -> str:
    """Create ColabFold-usable MSA for the monomer/multimer observation."""
    if '_' in symbol:  # multimer
        protein_a, protein_b = symbol.split('_')
        msa_a, msa_b = find_msa(protein_a, msa_dir), find_msa(protein_b, msa_dir)
        sequences_a = extract_header_sequence_pairs(f'{msa_dir}/{msa_a}')
        sequences_b = extract_header_sequence_pairs(f'{msa_dir}/{msa_b}')
        combined_msa_path = f'{msa_dir}/multimer/{symbol}.msa.a3m'
        write_combined_a3m(sequences_a, sequences_b, combined_msa_path)
        return combined_msa_path
    else:  # monomer
        msa = find_msa(symbol, msa_dir)
        return f'{msa_dir}/{msa}'


def create_observations(filepath, lower, upper) -> pd.DataFrame:
    """Create a dataframe of all required observations."""
    dataframe = pd.read_csv(filepath, names=['symbol', 'length'], skiprows=lower+1, nrows=upper-lower+1)
    logging.debug(dataframe.head(5))
    return dataframe


def run_colabfold_script(dataset_path, lower, upper, msa_dir) -> None:
    """Run the colabfold script for each observation."""
    df = create_observations(dataset_path, lower, upper)
    df['file'] = df['symbol'].apply(lambda symbol: prep_msas(symbol, msa_dir))
    logging.debug(df.head(5))
    

def parse_command_line():  # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-b',
        '--lower_bound',
        type=int,
        default=0,
        help='Which observation the program should start running ColabFold on (default: 0)'
    )
    parser.add_argument(
        '-u',
        '--upper_bound',
        type=int,
        default=10000,
        help='Which observation the program should finish running ColabFold on (default: 10,000)'
    )
    parser.add_argument(
        '-d',
        '--dataset_path',
        type=str,
        default='data/interim/sequence_lengths.csv',
        help='The path to the sorted dataset'
    )
    parser.add_argument(
        '-m',
        '--msa_dir',
        type=str,
        default='/cluster/projects/kumargroup/jackson/msas/',
        help='The path to the MSA files'
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
        else os.path.join(os.getcwd(), "logs/run_colabfold.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    if not os.path.exists(args.dataset_path):
        logging.warning(f'The dataset path {args.dataset_path} does not exist. Cancelling operation...')
        return
    if not os.path.exists(args.dataset_path):
        logging.warning(f'The MSA directory {args.msa_dir} does not exist. Cancelling operation...')
        return
    logging.info('Running colabfold script...')
    run_colabfold_script(args.dataset_path, args.lower_bound, args.upper_bound, args.msa_dir)

if __name__ == '__main__':  # pragma: no cover
    main()
