"""
A collection of functions to remove biased or erroneous
data observations from a dataset.
"""


import csv
import logging
import os

import pandas as pd


def chunk_input_genes(input_genes: list, chunk_size: int = 20) -> list:
    """Chunk input genes since the Biogrid API is limited to returning
    10,000 interactions."""
    chunked_list = [input_genes[i:i + chunk_size] for i in range(0, len(input_genes), chunk_size)]
    return chunked_list


def parse_input_genes(infile) -> list:
    """Parse the input file and return a list of official gene symbols."""
    df = pd.read_csv(infile)
    input_genes = df.iloc[:, 0].tolist()
    return input_genes


def remove_ground_truth_data(
        genes: list, 
        reference_file_path: str, 
        sheet_name: str,
        column_name: str
) -> list:
    """
    Remove genes from a dataset that were experimentally tested in 
    some reference dataset specified by an excel file.

    Parameters
    ----------
    unpruned_ppis: list 
        The list of PPIs, where each of form [gene_1]-[gene_2]
    reference_sheet_file_path: str
        The path to the reference sheet to be used to prune
        PPIs
    sheet_name: str
        The name of the sheet with gene names
    column_name: str
        The header for the column with gene names
    """
    # Read in all reference genes
    ground_truth_genes_df = pd.read_excel(
        reference_file_path, 
        sheet_name=sheet_name,
        usecols=[column_name],
        engine='calamine'
    )
    ground_truth_gene_set = set(ground_truth_genes_df[column_name])
    return [gene for gene in genes if gene not in ground_truth_gene_set]


def write_ppi_file(ppi_list, outfile):
    """Write the protein-protein interactions to a csv file with columns
    gene_name_a, gene_name_b where a is the gene of interest (e.g., cancer
    driver genes) and b is the interacting protein."""
    logging.debug(f'Writing to directory {os.path.abspath(outfile)}...')
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(outfile))), exist_ok=True)
    ppi_list.sort()
    rows = [pair.split('*') for pair in ppi_list]
    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['gene_symbol_a', 'gene_symbol_b'])
        writer.writerows(rows)
