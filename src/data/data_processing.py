"""
A collection of functions to remove biased or erroneous
data observations from a dataset.
"""


import pandas as pd


def remove_ground_truth_data(
        genes: list, 
        reference_file_path: str, 
        sheet_name: str,
        column_name: str,
        triplet_file: str
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
    triplet_file: str
        The path to the file of experimentally validated isoform-
        isoform interactions. Expected headers are ref_ID and bait_ID.
        Expects a CSV.
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


def parse_input_genes(infile) -> list:
    """Parse the input file and return a list of official gene symbols."""
    df = pd.read_csv(infile)
    input_genes = df.iloc[:, 0].tolist()
    return input_genes


def chunk_input_genes(input_genes: list, chunk_size: int = 20) -> list:
    """Chunk input genes since the Biogrid API is limited to returning
    10,000 interactions."""
    chunked_list = [input_genes[i:i + chunk_size] for i in range(0, len(input_genes), chunk_size)]
    return chunked_list
