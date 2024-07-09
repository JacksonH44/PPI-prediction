"""
A collection of functions that call bio database APIs and
process the resulting data.
"""

import logging
import os
import requests
import sys

import aiohttp

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from core import config as cfg


async def get_sequence_lengths(
    session: aiohttp.ClientSession, ensembl_transcript_ids: list[str]
) -> dict[str, int]:
    """Query the Ensembl API to get the number of amino acids in a sequence"""
    lookup = "/lookup/id"
    headers = {"Content-Type": "application/json"}
    data = {"ids": ensembl_transcript_ids, "expand": 1}
    logging.debug(f"Calling Ensembl API with {len(ensembl_transcript_ids)} transcripts")
    # response = requests.post(cfg.ENSEMBL_BASE_URL + lookup, headers=headers, json=data)
    request_url = cfg.ENSEMBL_BASE_URL + lookup
    response = await session.request(
        "POST", url=request_url, headers=headers, json=data
    )

    if response.status != 200:
        logging.warning(
            f"Unable to retrieve sequence for one of the following:\n{ensembl_transcript_ids}"
        )
        return {}

    res = await response.json()
    aa_counts = {}
    for transcript in ensembl_transcript_ids:
        if transcript in res and "Translation" in res[transcript]:
            aa_counts[transcript] = res[transcript]["Translation"]["length"]
    logging.debug(f"Found counts for {len(aa_counts.keys())} proteins...")
    return aa_counts


def find_uniprot_ids(transcripts: list[str]) -> dict[str, str]:
    """Query the Biomart API to map canonical Ensembl transcripts to
    UniProt IDs."""
    transcripts = [transcript.split(".")[0] for transcript in transcripts]
    transcript_ids_str = ",".join(transcripts)

    # The BioMart API requires a query in XML format
    # For parameter list go to:
    # https://www.ensembl.org/biomart/martservice?type=attributes&dataset=hsapiens_gene_ensembl
    query_xml = f"""
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="">
        <Dataset name="hsapiens_gene_ensembl" interface="default">
            <Filter name="ensembl_transcript_id" value="{transcript_ids_str}"/>
            <Attribute name="ensembl_transcript_id"/>
            <Attribute name="uniprotswissprot"/>
        </Dataset>
    </Query>
    """

    response = requests.post(
        "http://www.ensembl.org/biomart/martservice", data={"query": query_xml}
    )
    transcript_to_uniprot = {}
    if response.status_code == 200:
        lines = response.text.strip().split("\n")

        for line in lines[1:]:
            # Split TSV response and extract transcript ID
            # and UniProt SwissProt ID
            values = line.split("\t")
            transcript_id = values[0]
            try:
                uniprotswissprot_id = values[1]
                if uniprotswissprot_id == "":
                    logging.warning(
                        f"Was not able to find UniProt ID for protein with transcript ID {transcript_id}"
                    )
                else:
                    transcript_to_uniprot[transcript_id] = uniprotswissprot_id
            except IndexError:
                logging.warning(
                    f"Was not able to find UniProt ID for protein with transcript ID {transcript_id}"
                )
    else:
        logging.warning(
            f"Error: {response.status_code} for one of the transcripts in list:\n{transcripts}"
        )
    return transcript_to_uniprot


def filter_for_uniref30(proteins: list[str]) -> list[str]:
    """Filter out genes that don't appear in the Uniref30 DB from a list of PPIs or genes."""
    treat_as_ppis = False
    if "*" in proteins[0]:
        treat_as_ppis = True
        logging.debug("Uniref30 filter treating input as PPIs...")
        gene_set = set()
        for ppi in proteins:
            protein_a, protein_b = ppi.split("*")
            gene_set.add(protein_a)
            gene_set.add(protein_b)
            genes = list(gene_set)
    else:
        logging.debug("UniRef30 filter treating input as gene list...")
        genes = proteins
    gene_ids_str = ",".join(genes)

    # The BioMart API requires a query in XML format
    # For parameter list go to:
    # https://www.ensembl.org/biomart/martservice?type=attributes&dataset=hsapiens_gene_ensembl
    query_xml = f"""
    <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="">
        <Dataset name="hsapiens_gene_ensembl" interface="default">
            <Filter name="uniprot_gn_symbol" value="{gene_ids_str}"/>
            <Attribute name="uniprot_gn_symbol"/>
            <Attribute name="uniprotswissprot"/>
        </Dataset>
    </Query>
    """

    filtered_gene_set = set()
    response = requests.post(cfg.BIOMART_BASE_URL, data={"query": query_xml})
    if response.status_code == 200:
        lines = response.text.strip().split("\n")

        for line in lines[1:]:
            # Split TSV response and extract gene symbol ID
            # and UniProt SwissProt ID
            values = line.split("\t")
            gene_id = values[0]
            try:
                uniprotswissprot_id = values[1]
                if uniprotswissprot_id != "":
                    filtered_gene_set.add(gene_id)
            except IndexError:
                logging.warning(
                    f"""No SwissProt value for gene {gene_id}.
                    This could be an error in the API call, but ensure a SwissProt ID exists for this gene."""
                )
    filtered_genes = []
    if treat_as_ppis:
        for ppi in proteins:
            protein_a, protein_b = ppi.split("*")
            if protein_a in filtered_gene_set and protein_b in filtered_gene_set:
                filtered_genes.append(ppi)
    else:
        filtered_genes = list(filtered_gene_set)
    return filtered_genes


async def get_interactors(
    session: aiohttp.ClientSession,
    gene_list: list,
    cross_study_level: int = 2,
    relax_evidence: bool = False,
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
    relax_evidence: bool
        Whether to relax the criteria for a protein-protein interaction.
        If False, only high-confidence (low throughput, high reliability
        physical cross-referenced) evidence will be used. If True, any \
        evidence of a physical interaction in Biogrid will be used.
    """
    request_url = cfg.BIOGRID_BASE_URL + "/interactions"
    strict_evidence_list = cfg.BIOGRID_STRICT_EVIDENCE
    relaxed_evidence_list = cfg.BIOGRID_RELAXED_EVIDENCE
    evidence_list = (
        strict_evidence_list if not relax_evidence else relaxed_evidence_list
    )
    inter_species_excluded = "true" if not relax_evidence else "false"
    throughput_level = "low" if not relax_evidence else "any"
    params = {
        "accesskey": cfg.BIOGRID_API_KEY,
        "format": "json",
        "geneList": "|".join(gene_list),
        "interSpeciesExcluded": inter_species_excluded,
        "selfInteractionsExcluded": "true",
        "evidenceList": "|".join(evidence_list),
        "includeEvidence": "true",
        "includeInteractors": "true",
        "includeInteractorInteractions": "false",
        "searchNames": "true",
        "throughputTag": throughput_level,
        "taxId": 9606,  # Homo Sapiens taxonomy ID
        "includeHeader": "true",
    }
    resp = await session.request("GET", url=request_url, params=params)
    if resp.status != 200:
        logging.warning("Failed to get request for one of the genes in the gene list")
        return []

    interactions = await resp.json()
    if len(interactions) == 0:
        logging.warning(
            "Failed to retrieve any interaction experiments for the genes in the gene list"
        )
        return []

    # Create a hashmap of PPI pairs and the number of unique experiments. Unique experiments
    # are experiments that are from different studies (different Pubmed IDs).
    ppis: dict[str, set[str]] = {}
    for interaction in interactions.values():
        sym_a = interaction["OFFICIAL_SYMBOL_A"]
        sym_b = interaction["OFFICIAL_SYMBOL_B"]
        # Ensure pair is not already in the dataset
        experimental_id = (
            f'{interaction["PUBMED_ID"]}_{interaction["EXPERIMENTAL_SYSTEM"]}'
        )
        if f"{sym_a}*{sym_b}" in ppis:
            ppis[f"{sym_a}*{sym_b}"].add(experimental_id)
        elif f"{sym_b}*{sym_a}" in ppis:
            ppis[f"{sym_b}*{sym_a}"].add(experimental_id)
        # Write cancer driver gene first
        elif sym_a in gene_list:
            ppis[f"{sym_a}*{sym_b}"] = {experimental_id}
        else:  # Gene b in gene_list
            ppis[f"{sym_b}*{sym_a}"] = {experimental_id}
    # In the case of strict evidence, filter PPIs by number of cross-references
    if not relax_evidence:
        output_ppis = [
            ppi for ppi, e_id in ppis.items() if len(e_id) >= cross_study_level
        ]
    else:
        output_ppis = list(ppis.keys())
    logging.debug(
        f"Found a total of {len(interactions)} interactions and generated {len(output_ppis)} PPIs for:\n{gene_list}"
    )
    return output_ppis
