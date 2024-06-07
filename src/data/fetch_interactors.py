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


def get_interactors(gene_list):
    """
    Retrive a dataset of all interactors for all genes
    of interest.

    Parameters
    ----------
    gene_list: list(str)
        A list of the official symbols of a gene
    """
    request_url = cfg.BIOGRID_BASE_URL + "/interactions"
    evidence_list = [
        "affinity capture-luminescence",
        "affinity capture-ms",
        "affinity capture-rna",
        "affinity capture-western",
        "co-crystal structure",
        "two-hybrid"
    ]
    params = {
        "accesskey": cfg.BIOGRID_API_KEY,
        "format": "json",
        "geneList": "|".join(gene_list),
        "evidenceList": "|".join(evidence_list),
        "includeEvidence": "true",
        "searchNames": "true",
        "throughputTag": "low",
        "taxId": 9606, # Homo Sapiens taxonomy ID
        "includeHeader": "true"
    }

    r = requests.get(request_url, params=params)
    if r.status_code == 200:
        interactor_list = r.json()
        logging.debug(json.dumps(interactor_list, indent=2))
    else:
        logging.warning(f"Failed to get request for one of the genes in the gene list")


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
    get_interactors(["FANCE"])


if __name__ == '__main__':
    main()