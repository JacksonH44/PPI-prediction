"""
Fetch protein interactors for a set of genes.
"""


import argparse
import json
import logging
import os
import requests
import sys

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from core import config as cfg


def get_interactors(gene_list, cross_study_level=2):
    """
    Retrive a dataset of all interactors for all genes
    of interest.

    Parameters
    ----------
    gene_list: list(str)
        A list of the official symbols of a gene
    cross_study_level: int
        The minimum number of different studies that have 
        confirmed a PPI. Default is 2, and it is recommended
        that this value not be set to below 2
    """
    request_url = cfg.BIOGRID_BASE_URL + "/interactions"
    evidence_list = [
        "affinity capture-ms",
        "affinity capture-western",
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
    r = requests.get(request_url, params=params)
    if r.status_code != 200:
        logging.warning("Failed to get request for one of the genes in the gene list")
        return
    
    interactions = r.json()
    if len(interactions) == 0:
        logging.warning("Failed to retrieve any interaction experiments for the genes in the gene list")
        return
    
    logging.debug(f'Found a total of {len(interactions)} interactions...')
    logging.debug(json.dumps(interactions, indent=2))
    # Create a hashmap of PPI pairs and the number of unique experiments. Unique experiments
    # are experiments that are from different studies (different Pubmed IDs).
    ppis = {}
    for interaction in interactions.values():
        sym_a = interaction['OFFICIAL_SYMBOL_A']
        sym_b = interaction['OFFICIAL_SYMBOL_B']
        # Write cancer driver gene first
        pair = f'{sym_a}-{sym_b}' if sym_a in gene_list else f'{sym_b}-{sym_a}'
        if pair not in ppis:
            ppis[pair] = {interaction['PUBMED_ID']}
        else:
            ppis[pair].add(interaction['PUBMED_ID'])
    logging.debug(f'All PPIs:\n{ppis}')
    high_conf_ppis = [k for k, v in ppis.items() if len(v) >= cross_study_level]
    logging.debug(f'High confidence PPIs:\n{high_conf_ppis}')
    return high_conf_ppis


def parse_input_genes(infile):
    """Parse the input file and return a list of official gene symbols."""
    df = pd.read_csv(infile)
    input_genes = df.iloc[:, 0].tolist()
    logging.debug(input_genes)
    return input_genes


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
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
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/fetch_interactors.log")
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    parse_input_genes('data/raw/cancer_driver_gene_list.csv')
    # ppis = get_interactors(["TMPRSS2"], 2)


if __name__ == '__main__':
    main()