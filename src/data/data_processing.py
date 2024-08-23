"""
A collection of functions to remove biased or erroneous
data observations from a dataset.
"""

import csv
import logging
import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join("..", "..")))
from core import config as cfg


class UndersamplingError(Exception):
    """A class that represents an error in the undersampling process."""

    pass


def map_symbols_to_transcripts(
    symbols: list[str], reference_file: str = cfg.MANE_FILE
) -> dict[str, str]:
    """Map NCBI gene symbol to canonical transcript"""
    transcript_df = pd.read_csv(
        reference_file, sep="\t", usecols=["symbol", "Ensembl_nuc"]
    )
    filtered_df = transcript_df[transcript_df["symbol"].isin(symbols)]
    filtered_df.loc[:, "Ensembl_nuc"] = filtered_df["Ensembl_nuc"].apply(
        lambda x: x.split(".")[0]
    )
    return (filtered_df.set_index("symbol").to_dict())["Ensembl_nuc"]


def find_unique_genes(dataset_files) -> set[str]:
    """Find all unique genes in both the positive and negative datasets."""
    unique_genes: set[str] = set()
    for file_path in dataset_files:
        logging.debug(f"Processing file: {file_path}")
        with open(file_path, "r") as pos:
            pos_file = csv.reader(pos)
            next(pos_file)  # Skip header
            for row in pos_file:
                protein_a, protein_b = row[0], row[1]
                unique_genes.add(protein_a)
                unique_genes.add(protein_b)
    return unique_genes


def chunk_input_genes(input_genes: list, chunk_size: int = 20) -> list:
    """Chunk input genes since the Biogrid API is limited to returning
    10,000 interactions."""
    chunked_list = [
        input_genes[i : i + chunk_size] for i in range(0, len(input_genes), chunk_size)
    ]
    return chunked_list


def parse_input_genes(infile) -> list:
    """Parse the input file and return a list of official gene symbols."""
    df = pd.read_csv(infile)
    input_genes = df.iloc[:, 0].tolist()
    return input_genes


def remove_ground_truth_data(
    genes: list, reference_file_path: str, sheet_name: str, column_name: str
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
        engine="openpyxl",
    )
    ground_truth_gene_set = set(ground_truth_genes_df[column_name])
    return [gene for gene in genes if gene not in ground_truth_gene_set]


def write_ppi_file(ppi_list, outfile):
    """Write the protein-protein interactions to a csv file with columns
    gene_name_a, gene_name_b where a is the gene of interest (e.g., cancer
    driver genes) and b is the interacting protein."""
    logging.debug(f"Writing to directory {os.path.abspath(outfile)}...")
    os.makedirs(os.path.abspath(os.path.join(os.path.dirname(outfile))), exist_ok=True)
    ppi_list.sort()
    rows = [pair.split("*") for pair in ppi_list]
    with open(outfile, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["gene_symbol_a", "gene_symbol_b"])
        writer.writerows(rows)
