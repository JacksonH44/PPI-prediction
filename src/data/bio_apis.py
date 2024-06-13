"""
A collection of functions that call bio database APIs and 
process the resulting data.
"""


import logging
import os
import sys

import aiohttp

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../")))
from core import config as cfg


async def get_interactors(
    session: aiohttp.ClientSession,
    gene_list: list, 
    cross_study_level: int = 2,
    relax_evidence: bool = False
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
    evidence_list = strict_evidence_list if not relax_evidence else relaxed_evidence_list
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
        "taxId": 9606, # Homo Sapiens taxonomy ID
        "includeHeader": "true"
    }
    resp = await session.request('GET', url=request_url, params=params)
    if resp.status != 200:
        logging.warning("Failed to get request for one of the genes in the gene list")
        return
    
    interactions = await resp.json()
    if len(interactions) == 0:
        logging.warning("Failed to retrieve any interaction experiments for the genes in the gene list")
        return
    
    # Create a hashmap of PPI pairs and the number of unique experiments. Unique experiments
    # are experiments that are from different studies (different Pubmed IDs).
    ppis = {}
    for interaction in interactions.values():
        sym_a = interaction['OFFICIAL_SYMBOL_A']
        sym_b = interaction['OFFICIAL_SYMBOL_B']
        # Ensure pair is not already in the dataset
        experimental_id = f'{interaction["PUBMED_ID"]}_{interaction["EXPERIMENTAL_SYSTEM"]}'
        if f'{sym_a}*{sym_b}' in ppis:
            ppis[f'{sym_a}*{sym_b}'].add(experimental_id)
        elif f'{sym_b}*{sym_a}' in ppis:
            ppis[f'{sym_b}*{sym_a}'].add(experimental_id)
        # Write cancer driver gene first
        elif sym_a in gene_list:
            ppis[f'{sym_a}*{sym_b}'] = {experimental_id}
        else: # Gene b in gene_list
            ppis[f'{sym_b}*{sym_a}'] = {experimental_id}
    # In the case of strict evidence, filter PPIs by number of cross-references
    if not relax_evidence:
        output_ppis = [ppi for ppi, e_id in ppis.items() if len(e_id) >= cross_study_level]
    else:
        output_ppis = ppis.keys()
    logging.debug(f'Found a total of {len(interactions)} interactions and generated {len(output_ppis)} PPIs for:\n{gene_list}')
    return output_ppis