"""
Generate negative dataset for PPIs.
"""


import argparse
import asyncio
from collections import defaultdict
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


def find_subcellular_proteins(locations_df) -> defaultdict:
    """Find all potential proteins that are in the same
    subcellular location as a gene."""
    locations_to_genes = defaultdict(set)
    for _, row in locations_df.iterrows():
        gene = row['Gene name']
        locations = row['Main location']
        for location in locations.split(';'):
            locations_to_genes[location].add(gene)
    
    genes_sharing_locations = defaultdict(set)
    for genes in locations_to_genes.values():
        for gene in genes:
            # Add the gene itself to the off limits proteins to not allow for self-interactions
            genes_sharing_locations[gene] = genes_sharing_locations[gene] | genes
    return genes_sharing_locations


def find_interacting_proteins(all_ppis) -> defaultdict:
    """Find all proteins that have been recorded to 
    interact with a protein produced by the gene."""
    protein_dict = defaultdict(set)
    for pair in all_ppis:
        protein_a, protein_b = pair.split('*')
        protein_dict[protein_a].add(protein_b)
        protein_dict[protein_b].add(protein_a)
    return protein_dict


def find_unsuitable_partners(
        genes: list,
        locations_df: pd.DataFrame,
        all_ppis: list
) -> defaultdict:
    """
    Find all unsuitable partners for a list of genes and return
    a mapping of a gene to its unsuitable partners. By finding 
    unsuitable protein partners for a gene we save space as the 
    number of unsuitable partners for a gene is likely much less than
    the number of suitable partners.
    
    Parameters
    ----------
    genes : list
        A list of the genes of interest
    locations_df : pd.DataFrame
        A dataframe of all genes you are considering for the negative set
        and their subcellular locations
    all_ppis: list
        A list of all recorded protein-protein interactions involving at least
        one of the genes in the genes of interest list
    
    Returns
    -------
    unsuitable_partners : defaultdict
        key-value pairs where the gene of interest is the key and all the 
        proteins that are unsuitable to use as negative cases in a dataset 
        (in the same subcellular location or involved in a recorded PPI with
        the gene of interest).
    """
    unsuitable_partners = defaultdict(set)
    interacting_proteins = find_interacting_proteins(all_ppis)
    logging.debug(f'Found interacting proteins for {len(interacting_proteins.keys())} genes...')
    same_subcellular_proteins = find_subcellular_proteins(locations_df)
    logging.debug(f'''Found proteins with the same subcellular location for 
                  {len(same_subcellular_proteins.keys())} genes...''')
    for gene in genes:
        unsuitable_partners[gene] = interacting_proteins[gene] | same_subcellular_proteins[gene] # Set union
    return unsuitable_partners


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
    # Select for only reliablly found gene locations and have a main location
    location_df = location_df[(location_df['Reliability'] != 'Uncertain') & (location_df['Main location'].notna())]
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
    logging.debug('Finding unsuitable partners for each gene...')
    unsuitable_partners = find_unsuitable_partners(ground_truth_out_genes, locations_df, all_ppis)
    logging.debug(f'Found unsuitable partners for {len(unsuitable_partners)} genes...')


if __name__ == '__main__': # pragma: no cover
    asyncio.run(main())