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


def write_ppi_file(ppi_list, outfile):
    """Write the protein-protein interactions to a csv file with columns
    gene_name_a, gene_name_b where a is the gene of interest (e.g., cancer
    driver genes) and b is the interacting protein."""
    logging.debug(f'Writing to directory {os.path.abspath(outfile)}...')
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(outfile))), exist_ok=True)
    rows = [pair.split('-') for pair in ppi_list]
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
        "reconstituted-complex"
        "two-hybrid",
        "pca"
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
    # r = requests.get(request_url, params=params)
    resp = await session.request('GET', url=request_url, params=params)
    if resp.status != 200:
        logging.warning("Failed to get request for one of the genes in the gene list")
        return
    
    interactions = await resp.json()
    if len(interactions) == 0:
        logging.warning("Failed to retrieve any interaction experiments for the genes in the gene list")
        return
    
    logging.debug(f'Found a total of {len(interactions)} interactions...')
    # Create a hashmap of PPI pairs and the number of unique experiments. Unique experiments
    # are experiments that are from different studies (different Pubmed IDs).
    ppis = {}
    for interaction in interactions.values():
        sym_a = interaction['OFFICIAL_SYMBOL_A']
        sym_b = interaction['OFFICIAL_SYMBOL_B']
        # Ensure pair is not already in the dataset
        if f'{sym_a}-{sym_b}' in ppis:
            ppis[f'{sym_a}-{sym_b}'].add(interaction['PUBMED_ID'])
        elif f'{sym_b}-{sym_a}' in ppis:
            ppis[f'{sym_b}-{sym_a}'].add(interaction['PUBMED_ID'])
        # Write cancer driver gene first
        elif sym_a in gene_list:
            ppis[f'{sym_a}-{sym_b}'] = {interaction['PUBMED_ID']}
        else: # Gene b in gene_list
            ppis[f'{sym_b}-{sym_a}'] = {interaction['PUBMED_ID']}
    high_conf_ppis = [k for k, v in ppis.items() if len(v) >= cross_study_level]
    logging.debug(f'High confidence PPIs:\n{high_conf_ppis}')
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
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/fetch_interactors.log")
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
        tasks = []
        for chunk in chunked_genes:
            tasks.append(get_interactors(session, chunk, 2))
        ppi_lists = await asyncio.gather(*tasks, return_exceptions=True)
    finish = time.perf_counter()
    logging.debug(f'Finished curation in {round(finish - start, 2)} second(s)')
    ppis = [ppi for sublist in ppi_lists for ppi in sublist]
    logging.debug(f'Curated a total of {len(ppis)} PPIs')
    write_ppi_file(ppis, args.outfile)


if __name__ == '__main__': # pragma: no cover
    asyncio.run(main())