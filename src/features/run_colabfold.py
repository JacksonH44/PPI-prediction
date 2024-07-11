"""
Run ColabFold for protein pairs in a dataset.
"""

import argparse
import logging
import os

import pandas as pd


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


def create_observations(filepath, lower, upper) -> pd.DataFrame:
    """Create a dataframe of all required observations."""
    dataframe = pd.read_csv(filepath, names=['gene_symbol_a', 'gene_symbol_b'], skiprows=lower+1, nrows=upper-lower+1)
    dataframe['file'] = dataframe['gene_symbol_a'] + '_' + dataframe['gene_symbol_b']
    return dataframe


def run_colabfold_script(lower, upper, pos_path, neg_path, msa_dir) -> None:
    """Run the colabfold script for each observation."""
    pos_dataset = create_observations(pos_path, lower, upper)
    neg_dataset = create_observations(neg_path, lower, upper)
    logging.debug('Finding MSA files for all observations...')
    pos_dataset['gene_a_file'] = pos_dataset['gene_symbol_a'].apply(lambda gene_symbol: find_msa(gene_symbol, msa_dir))
    pos_dataset['gene_b_file'] = pos_dataset['gene_symbol_b'].apply(lambda gene_symbol: find_msa(gene_symbol, msa_dir))
    neg_dataset['gene_a_file'] = neg_dataset['gene_symbol_a'].apply(lambda gene_symbol: find_msa(gene_symbol, msa_dir))
    neg_dataset['gene_b_file'] = neg_dataset['gene_symbol_b'].apply(lambda gene_symbol: find_msa(gene_symbol, msa_dir))
    logging.debug(pos_dataset.head(3))
    logging.debug(neg_dataset.head(3))


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
        '-p',
        '--positive_dataset',
        type=str,
        default='data/processed/positive_ppis.csv',
        help='The path to the positive dataset'
    )
    parser.add_argument(
        '-n',
        '--negative_dataset',
        type=str,
        default='data/processed/negative_ppis.csv',
        help='The path to the negative dataset'
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
    run_colabfold_script(args.lower_bound, args.upper_bound, args.positive_dataset, args.negative_dataset, args.msa_dir)

if __name__ == '__main__':  # pragma: no cover
    main()
