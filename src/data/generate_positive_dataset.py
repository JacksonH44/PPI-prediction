"""
Fetch protein interactors for a set of genes.
"""


import argparse
import asyncio
import csv
import logging
import os
import sys
import time

import pandas as pd
import aiohttp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from core import config as cfg
from src.data.data_pruning import remove_ground_truth_data


def write_ppi_file(ppi_list, outfile):
    """Write the protein-protein interactions to a csv file with columns
    gene_name_a, gene_name_b where a is the gene of interest (e.g., cancer
    driver genes) and b is the interacting protein."""
    logging.debug(f'Writing to directory {os.path.abspath(outfile)}...')
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(outfile))), exist_ok=True)
    rows = [pair.split('_') for pair in ppi_list]
    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['gene_symbol_a', 'gene_symbol_b'])
        writer.writerows(rows)


async def get_interactors(
    session: aiohttp.ClientSession,
    gene_list: list, 
    cross_study_level: int = 2
) -> list:
    """
    Retrive a dataset of all interactors for all genes
    of interest.

    Parameters
    ----------
    session: aiohttp.ClientSession
        The current asynchronous session to make the async 
        calls to the BioGrid API.
    gene_list: list(str)
        A list of the official symbols of a gene
    cross_study_level: int
        The minimum number of different studies that have 
        confirmed a PPI. Default is 2, and it is recommended
        that this value not be set to below 2
    """
    request_url = cfg.BIOGRID_BASE_URL + "/interactions"
    evidence_list = [
        "affinity capture-luminescence",
        "affinity capture-ms",
        "affinity capture-rna",
        "affinity capture-western",
        "two-hybrid"
    ]
    params = {
        "accesskey": cfg.BIOGRID_API_KEY,
        "format": "json",
        "geneList": "|".join(gene_list),
        "interSpeciesExcluded": "true",
        "selfInteractionsExcluded": "true",
        "evidenceList": "|".join(evidence_list),
        "includeEvidence": "true",
        "includeInteractors": "true",
        "searchNames": "true",
        "throughputTag": "low",
        "taxId": 9606, # Homo Sapiens taxonomy ID
        "includeHeader": "true"
    }
    resp = await session.request('GET', url=request_url, params=params)
    if resp.status != 200:
        logging.warning("Failed to get request for one of the genes in the gene list")
        return
    
    interactions = await resp.json()
    if len(interactions) == 0:
        logging.warning("Failed to retrieve any interaction experiments for the genes in the gene list")
        return
    
    # Create a hashmap of PPI pairs and the number of unique experiments. Unique experiments
    # are experiments that are from different studies (different Pubmed IDs).
    ppis = {}
    for interaction in interactions.values():
        sym_a = interaction['OFFICIAL_SYMBOL_A']
        sym_b = interaction['OFFICIAL_SYMBOL_B']
        # Ensure pair is not already in the dataset
        experimental_id = f'{interaction["PUBMED_ID"]}_{interaction["EXPERIMENTAL_SYSTEM"]}'
        if f'{sym_a}_{sym_b}' in ppis:
            ppis[f'{sym_a}_{sym_b}'].add(experimental_id)
        elif f'{sym_b}_{sym_a}' in ppis:
            ppis[f'{sym_b}_{sym_a}'].add(experimental_id)
        # Write cancer driver gene first
        elif sym_a in gene_list:
            ppis[f'{sym_a}_{sym_b}'] = {experimental_id}
        else: # Gene b in gene_list
            ppis[f'{sym_b}_{sym_a}'] = {experimental_id}
    high_conf_ppis = [ppi for ppi, e_id in ppis.items() if len(e_id) >= cross_study_level]
    logging.debug(f'Found a total of {len(interactions)} and generated {len(high_conf_ppis)} high confidence PPIs...')
    return high_conf_ppis


def chunk_input_genes(input_genes: list, chunk_size: int = 60) -> list:
    """Chunk input genes since the Biogrid API is limited to returning
    10,000 interactions."""
    chunked_list = [input_genes[i:i + chunk_size] for i in range(0, len(input_genes), chunk_size)]
    return chunked_list


def parse_input_genes(infile) -> list:
    """Parse the input file and return a list of official gene symbols."""
    df = pd.read_csv(infile)
    input_genes = df.iloc[:, 0].tolist()
    logging.debug(input_genes)
    return input_genes


def parse_command_line(): # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Path to Input CSV containing genes of interest')
    parser.add_argument('outfile', type=str,
                        help='Output file path to write list of genes to')
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
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    input_genes = parse_input_genes(args.infile)
    chunked_genes = chunk_input_genes(input_genes)
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
    ground_truth_out_ppis = remove_ground_truth_data(
        ppis,
        cfg.GROUND_TRUTH_PATH,
        cfg.GROUND_TRUTH_SHEET,
        cfg.GROUND_TRUTH_COLUMN
    )
    logging.info(f'Curated a total of {len(ground_truth_out_ppis)} unbiased PPIs...')
    write_ppi_file(ground_truth_out_ppis, args.outfile)


if __name__ == '__main__': # pragma: no cover
    asyncio.run(main())