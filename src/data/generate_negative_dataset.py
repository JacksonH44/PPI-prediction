"""
Generate negative dataset for PPIs.
"""

import argparse
import asyncio
from collections import defaultdict
import csv
import logging
import os
import random
import sys
import time

import aiohttp
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from core import config as cfg
from src.data.bio_apis import get_interactors, get_sequence_lengths
from src.data.data_processing import (
    chunk_input_genes,
    map_symbols_to_transcripts,
    parse_input_genes,
    remove_ground_truth_data,
    UndersamplingError,
    write_ppi_file,
)

# Do random seeding
random.seed(cfg.SEED)


def randomly_select_partners(available_partners, num_samples):
    """Randomly select num_samples partners from the available partners."""
    if num_samples > len(available_partners):
        raise UndersamplingError
    partner_list = list(available_partners)
    sample = random.sample(partner_list, num_samples)
    return sample


def count_gene_symbols(positive_ppis):
    """Count the number of each gene symbol in the positve dataset."""
    gene_count = defaultdict(int)
    # Checking the file exists
    if not os.path.isfile(positive_ppis):
        raise FileNotFoundError(
            f"""The file {positive_ppis} does not exist.
                                Have you created the positive dataset first
                                (generate_positive_dataset.py)?"""
        )
    with open(positive_ppis, "r") as pos:
        csv_reader = csv.DictReader(pos)
        for row in csv_reader:
            gene_count[row["gene_symbol_a"]] += 1
    return dict(gene_count)


def undersample_dataset(
    locations_df: pd.DataFrame, positive_ppis: str, unsuitable_partners: dict
) -> list:
    """
    Create a list of negative PPIs where the number of PPIs for each
    gene of interest is the same as in the positive dataset.

    Parameters
    ----------
    locations_df : pd.DataFrame
        A dataframe of all potential partners a protein could have
    positive_ppis : str
        Path to the file holding all positive PPIs
    unsuitable_partners : dict
        A map of gene : unsuitable_partner pairs

    Returns
    -------
    neg_ppis : list
        A list of negative PPIs that has the same number of PPIs for
        each protein of interest as the positive set
    """
    neg_ppis = []
    all_proteins = set(locations_df["Gene name"])
    gene_count = count_gene_symbols(positive_ppis)
    for gene, num_ppis in gene_count.items():
        available_partners = all_proteins - unsuitable_partners[gene]
        # Check that there are enough available partners
        try:
            partners = randomly_select_partners(available_partners, num_ppis)
            ppis = [f"{gene}*{partner}" for partner in partners]
            neg_ppis.extend(ppis)
        except UndersamplingError:
            msg = f"""gene {gene} is trying to undersample with {num_ppis}
            PPIS but only has {len(available_partners)}. {gene} is not
            included in the negative dataset."""
            logging.warning(msg)
    return neg_ppis


async def find_too_long_proteins(
    genes: list[str], locations_df: pd.DataFrame
) -> defaultdict[str, set[str]]:
    """
    Find proteins that are too long to make complexes with each gene

    Parameters
    ----------
    genes : list[str]
        A list of all genes of interest
    locations_df : pd.DataFrame
        A dataframe of all possible genes that can be interacted with

    Returns
    -------
    length_map : defaultdict[str, set[str]]
        A hashmap of (gene name, proteins with too long sequence) pairs
    """
    length_map: defaultdict[str, set[str]] = defaultdict(set)
    all_genes = locations_df["Gene name"].to_list()
    all_transcripts = map_symbols_to_transcripts(all_genes)
    genes_of_interest_transcripts = map_symbols_to_transcripts(genes)
    logging.debug(
        f"Found {len(genes_of_interest_transcripts.keys())} mappings for genes of interest..."
    )
    all_transcripts = {**all_transcripts, **genes_of_interest_transcripts}
    logging.debug(
        f"Found {len(all_transcripts.keys())} genes to find lengths for\nHead:\n{list(all_transcripts.keys())[:5]}"
    )
    chunked_transcripts = chunk_input_genes(list(all_transcripts.values()), 225)
    async with aiohttp.ClientSession() as session:
        tasks = [get_sequence_lengths(session, chunk) for chunk in chunked_transcripts]
        mapped_sequences = await asyncio.gather(*tasks, return_exceptions=True)
    if isinstance(mapped_sequences, BaseException):
        logging.warning("Error in async API call. Exiting...")
        return length_map
    transcript_map = {
        transcript: length
        for hashmap in mapped_sequences
        for transcript, length in hashmap.items()  # type: ignore
    }
    for gene in genes:
        try:
            length_map[gene] = set(
                [
                    symbol
                    for symbol, transcript in all_transcripts.items()
                    if transcript_map[transcript]
                    + transcript_map[all_transcripts[gene]]
                    > 2000
                ]
            )
            logging.debug(f"Found {len(length_map[gene])} unworthy proteins for {gene}")
        except KeyError:
            logging.warning(f"Could not find transcript map for {gene}")

    return length_map


def find_subcellular_proteins(locations_df) -> defaultdict[str, set[str]]:
    """Find all potential proteins that are in the same
    subcellular location as a gene."""
    locations_to_genes = defaultdict(set)
    for _, row in locations_df.iterrows():
        gene = row["Gene name"]
        locations = row["Main location"]
        for location in locations.split(";"):
            locations_to_genes[location].add(gene)

    genes_sharing_locations: defaultdict[str, set[str]] = defaultdict(set)
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
        protein_a, protein_b = pair.split("*")
        protein_dict[protein_a].add(protein_b)
        protein_dict[protein_b].add(protein_a)
    return protein_dict


async def find_unsuitable_partners(
    genes: list, locations_df: pd.DataFrame, all_ppis: list
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
    logging.debug(
        f"Found interacting proteins for {len(interacting_proteins.keys())} genes..."
    )
    same_subcellular_proteins = find_subcellular_proteins(locations_df)
    logging.debug(
        f"""Found proteins with the same subcellular location for
                  {len(same_subcellular_proteins.keys())} genes..."""
    )
    logging.debug(
        "Finding proteins that make complexes that are too long for folding..."
    )
    too_long_proteins = await find_too_long_proteins(genes, locations_df)
    for gene in genes:
        unsuitable_partners[gene] = (
            interacting_proteins[gene]
            | same_subcellular_proteins[gene]
            | too_long_proteins[gene]
        )  # Set union
    return unsuitable_partners


def get_locations(gene_file: str, location_file: str) -> pd.DataFrame:
    """Get subcellular locations for all genes in the file."""
    gene_df = pd.read_csv(gene_file, sep="\t", usecols=["symbol"])
    gene_df = gene_df.rename(columns={"symbol": "Gene name"})
    # Rename column to be able to merge dataframes
    logging.debug(f"Number of genes initially: {gene_df.shape[0]}")
    location_df = pd.read_csv(
        location_file, sep="\t", usecols=["Gene name", "Reliability", "Main location"]
    )
    # Select for only reliablly found gene locations and have a main location
    location_df = location_df[
        (location_df["Reliability"] != "Uncertain")
        & (location_df["Main location"].notna())
    ]
    logging.debug(f"Head of location dataframe:\n{location_df.head()}")
    logging.debug("Merging dataframes...")
    merged_df = gene_df.merge(location_df, on="Gene name", how="inner", copy=False)
    return merged_df


def parse_command_line():  # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-g",
        "--gene_file",
        type=str,
        default="data/raw/cancer_driver_gene_list.csv",
        help="Path to CSV containing genes of interest",
    )
    parser.add_argument(
        "-a",
        "--all_genes",
        type=str,
        default="data/raw/MANE_summary_v3.csv",
        help="Path to CSV containing all genes (could be MANE file)",
    )
    parser.add_argument(
        "-c",
        "--location_file",
        type=str,
        default="data/raw/subcellular_location.csv",
        help="Path to TSV containing subcellular locations for all genes",
    )
    parser.add_argument(
        "-p",
        "--positive_dataset",
        type=str,
        default="data/processed/positive_ppis.csv",
        help="Path to CSV file fo positive PPI dataset",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="data/processed/negative_ppis.csv",
        help="Path to output file to store negative PPI datset",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Change logging level from default level to noisiest level",
    )
    parser.add_argument(
        "-l",
        "--logfile",
        type=str,
        default=None,
        help="Specify the name of the log file",
    )
    args = parser.parse_args()
    return args


async def main():  # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    logfile = (
        args.logfile
        if args.logfile is not None
        else os.path.join(os.getcwd(), "logs/generate_negative_dataset.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    input_genes = parse_input_genes(args.gene_file)
    logging.debug("Removing genes with no MANE transcript...")
    mane_df = pd.read_csv(cfg.MANE_FILE, sep="\t", usecols=["symbol"])
    input_genes = [gene for gene in input_genes if gene in mane_df["symbol"].values]
    ground_truth_out_genes = remove_ground_truth_data(
        input_genes,
        cfg.GROUND_TRUTH_PATH,
        cfg.GROUND_TRUTH_SHEET,
        cfg.GROUND_TRUTH_COLUMN,
    )
    logging.debug(
        f"Found {len(ground_truth_out_genes)} genes from {len(input_genes)} after removing ground truth genes..."
    )
    # Generate a list of all possible proteins a protein of interest could interact with
    # and get subcellular location
    logging.debug("Getting subcellular locations...")
    locations_df = get_locations(args.all_genes, args.location_file)
    logging.debug(
        f"Found {locations_df.shape[0]} potential protein partners for each gene (before filtering)..."
    )
    logging.debug("Generating all protein partners for each gene of interest...")
    chunked_genes = chunk_input_genes(ground_truth_out_genes, 10)
    start = time.perf_counter()  # Time the function call for debugging
    async with aiohttp.ClientSession() as session:
        tasks = [
            get_interactors(session, chunk, 2, relax_evidence=True)
            for chunk in chunked_genes
        ]
        ppi_lists = await asyncio.gather(*tasks, return_exceptions=True)
    all_ppis = [ppi for sublist in ppi_lists for ppi in sublist]
    finish = time.perf_counter()
    logging.info(f"Finished curation in {round(finish - start, 2)} second(s)")
    logging.info(f"Curated a total of {len(all_ppis)} PPIs...")
    logging.debug("Finding unsuitable partners for each gene...")
    unsuitable_partners = await find_unsuitable_partners(
        ground_truth_out_genes, locations_df, all_ppis
    )
    logging.debug(f"Found unsuitable partners for {len(unsuitable_partners)} genes...")
    logging.debug(
        f"""EXAMPLE: Unsuitable partners for HLF
                  ({len(unsuitable_partners["HLF"])} total)"""
    )
    logging.debug("Undersampling to create negative dataset...")
    try:
        neg_ppis = undersample_dataset(
            locations_df, args.positive_dataset, unsuitable_partners
        )
        logging.debug(f"Found {len(neg_ppis)} negative PPIs...")
        logging.debug(f"Writing negative PPIs to {args.outfile}...")
        write_ppi_file(neg_ppis, args.outfile)
    except FileNotFoundError:
        msg = f"The file {args.positive_dataset} does not exist. The negative dataset was not processed."
        logging.warning(msg)
    except UndersamplingError as e:
        logging.warning(e)


if __name__ == "__main__":  # pragma: no cover
    asyncio.run(main())
