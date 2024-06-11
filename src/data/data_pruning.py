"""A collection of functions to remove biased or erroneous
data observations from a dataset"""


import pandas as pd


def remove_ground_truth_data(
        unpruned_ppis: list, 
        reference_file_path: str, 
        sheet_name: str,
        column_name: str
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
    """
    def is_interaction_valid(interaction, gene_set) -> bool:
        """Determines whether a given PPI from the unpruned
        PPIs contains any proteins that were tested in the
        ground truth dataset."""
        protein_1, protein_2 = interaction.split('-')
        return protein_1 not in gene_set and protein_2 not in gene_set
    
    genes_df = pd.read_excel(
        reference_file_path, 
        sheet_name=sheet_name,
        usecols=[column_name]
    )
    gene_set = set(genes_df[column_name])
    filtered_interactions = [interaction for interaction in unpruned_ppis if is_interaction_valid(interaction, gene_set)]
    return filtered_interactions
