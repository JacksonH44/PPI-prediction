"""
Generate negative dataset for PPIs.
"""


import argparse
import asyncio
import logging
import os
import sys
import time

import aiohttp
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from core import config as cfg
from src.data.bio_apis import get_interactors
from src.data.data_processing import chunk_input_genes, parse_input_genes, remove_ground_truth_data


def get_locations(gene_file : str, location_file : str) -> pd.DataFrame:
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
    return merged_df


def parse_command_line(): # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene_file', type=str,
                        help='Path to CSV containing genes of interest')
    parser.add_argument('all_genes', type=str,
                        help='Path to CSV containing all genes (could be MANE file)')
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


async def main(): # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/generate_negative_dataset.log")
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    input_genes = parse_input_genes(args.gene_file)
    ground_truth_out_genes = remove_ground_truth_data(
        input_genes,
        cfg.GROUND_TRUTH_PATH,
        cfg.GROUND_TRUTH_SHEET,
        cfg.GROUND_TRUTH_COLUMN,
        cfg.TRIPLET_FILE
    )
    logging.debug(f'Found {len(ground_truth_out_genes)} genes from {len(input_genes)} after removing ground truth genes...')
    # Generate a list of all possible proteins a protein of interest could interact with
    # and get subcellular location
    logging.debug('Getting subcellular locations...')
    locations_df = get_locations(args.all_genes, args.location_file)
    logging.debug(f'Found {locations_df.shape[0]} potential protein partners for each gene (before filtering)...')
    logging.debug('Generating all protein partners for each gene of interest...')
    chunked_genes = chunk_input_genes(ground_truth_out_genes, 10)
    start = time.perf_counter() # Time the function call for debugging
    async with aiohttp.ClientSession() as session:
        tasks = [get_interactors(session, chunk, 2, relax_evidence=True) for chunk in chunked_genes]
        ppi_lists = await asyncio.gather(*tasks, return_exceptions=True)
    all_ppis = [ppi for sublist in ppi_lists for ppi in sublist]
    finish = time.perf_counter()
    logging.info(f'Finished curation in {round(finish - start, 2)} second(s)')
    logging.info(f'Curated a total of {len(all_ppis)} PPIs...')


if __name__ == '__main__': # pragma: no cover
    asyncio.run(main())