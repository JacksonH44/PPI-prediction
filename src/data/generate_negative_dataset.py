"""
Generate negative dataset for PPIs.
"""


import argparse
import logging
import os

import pandas as pd


def get_locations(gene_file, location_file):
    """Get subcellular locations for all genes in the file."""
    gene_df = pd.read_csv(
        gene_file,
        usecols=['gene']
    )
    gene_df.columns = ['Gene name']
    # Rename column to be able to merge dataframes
    logging.debug(f'Number of genes initially: {gene_df.shape[0]}')
    location_df = pd.read_csv(
        location_file,
        sep='\t',
        usecols=['Gene name', 'Reliability', 'Main location']
    )
    # Select for only reliablly found gene locations
    location_df = location_df[location_df['Reliability'] != 'Uncertain']
    logging.debug(f'Head of location dataframe:\n{location_df.head()}')
    logging.debug('Merging dataframes...')
    merged_df = gene_df.merge( 
        location_df, 
        on='Gene name',
        how='inner',
        copy=False
    )
    logging.debug(f'Merged dataframe:\n{merged_df.head()}')
    logging.debug(f'Merged dataframe is of size {merged_df.shape[0]}')


def parse_command_line(): # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene_file', type=str,
                        help='Path to CSV containing genes of interest')
    parser.add_argument('location_file', type=str,
                        help='Path to TSV containing subcellular locations for all genes')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Change logging level from default level to noisiest level')
    parser.add_argument('-l', '--logfile',
                        type=str, default=None,
                        help='Specify the name of the log file')
    args = parser.parse_args()
    return args


def main(): # pragma : no cover
    """Run the command line program."""
    args = parse_command_line()
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/generate_negative_dataset.log")
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    logging.debug("Getting subcellular locations...")
    get_locations(args.gene_file, args.location_file)


if __name__ == '__main__': # pragma : no cover
    main()