"""
A collection of functions to remove biased or erroneous
data observations from a dataset.
"""


import pandas as pd


def remove_ground_truth_data(
        unpruned_ppis: list, 
        reference_file_path: str, 
        sheet_name: str,
        column_name: str,
        triplet_file: str
) -> list:
    """
    Remove PPIs from a dataset that were tested in Y2H
    in the ground truth isoform dataset.

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
    def is_interaction_valid(interaction, gene_set, interaction_set) -> bool:
        """Determines whether a given PPI from the unpruned
        PPIs contains any proteins that were tested in the
        ground truth dataset."""
        protein_1, protein_2 = interaction.split('_')
        return protein_1 not in gene_set and protein_2 not in gene_set and interaction not in interaction_set
    
    # Read in all reference genes
    genes_df = pd.read_excel(
        reference_file_path, 
        sheet_name=sheet_name,
        usecols=[column_name],
        engine='calamine'
    )
    gene_set = set(genes_df[column_name])
    # Read in all experimentally validated interactions
    experimental_interactions = pd.read_csv(
        triplet_file,
        usecols = ['ref_ID', 'bait_ID']
    )
    ref_bait_strings = experimental_interactions['ref_ID'].astype(str) + '_' + experimental_interactions['bait_ID'].astype(str)
    bait_ref_strings = experimental_interactions['bait_ID'].astype(str) + '_' + experimental_interactions['ref_ID'].astype(str)
    interaction_set = set(ref_bait_strings).union(set(bait_ref_strings))
    filtered_interactions = [interaction for interaction in unpruned_ppis if is_interaction_valid(interaction, gene_set, interaction_set)]
    return filtered_interactions


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
