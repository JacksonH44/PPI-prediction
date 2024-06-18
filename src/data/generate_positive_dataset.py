"""
Generate positive dataset for PPIs.
"""


import argparse
import asyncio
import logging
import os
import time

import aiohttp

from core import config as cfg
from src.data.bio_apis import get_interactors
from src.data.data_processing import (
    chunk_input_genes,
    parse_input_genes,
    remove_ground_truth_data,
    write_ppi_file
)


def parse_command_line(): # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Path to Input CSV containing genes of interest')
    parser.add_argument('-o', '--outfile', type=str, 
                        default='data/processed/positive_ppis.csv',
                        help='Output file path to write positive dataset to')
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
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/generate_positive_dataset.log")
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    input_genes = parse_input_genes(args.infile)
    ground_truth_out_genes = remove_ground_truth_data(
        input_genes,
        cfg.GROUND_TRUTH_PATH,
        cfg.GROUND_TRUTH_SHEET,
        cfg.GROUND_TRUTH_COLUMN
    )
    logging.debug(f'Found {len(ground_truth_out_genes)} genes from {len(input_genes)} after removing ground truth genes...')
    chunked_genes = chunk_input_genes(ground_truth_out_genes)
    start = time.perf_counter() # Time the function call for debugging
    # Make the API calls asynchronously
    async with aiohttp.ClientSession() as session:
        tasks = [get_interactors(session, chunk, 2) for chunk in chunked_genes]
        ppi_lists = await asyncio.gather(*tasks, return_exceptions=True)
    finish = time.perf_counter()
    logging.info(f'Finished curation in {round(finish - start, 2)} second(s)')
    ppis = [ppi for sublist in ppi_lists for ppi in sublist]
    logging.info(f'Curated a total of {len(ppis)} high confidence PPIs...')
    logging.debug('Removing ground truth PPIs...')
    # Remove PPIs that include a protein from the ground truth data
    logging.info(f'Curated a total of {len(ppis)} unbiased PPIs...')
    write_ppi_file(ppis, args.outfile)


if __name__ == '__main__': # pragma: no cover
    asyncio.run(main())
