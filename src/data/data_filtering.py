"""
Filtering operations on datasets.
"""

import time
import logging

import aiohttp

from src.data.bio_apis import get_sequence_lengths
from src.data.data_processing import map_symbols_to_transcripts


async def filter_out_long_sequences(
    session: aiohttp.ClientSession, ppis: list[str]
) -> list[str]:
    """
    Filter out sequences that make a protein
    pair greater than 2000 amino acids.
    
    Parameters
    ----------
    session : aiohttp.ClientSession
        The Session object that holds the context in which this function is called
        asynchronously
    ppis : list[str]
        A list of protein-protein interactions (PPIs) from which long sequences
        need to be filtered out from
    
    Returns
    -------
    filtered_ppis : list[str]
        The input PPI list with all PPIs with total sequence greater than 2,000 amino
        acids filtered out
    """
    ppi_symbols = []
    for ppi in ppis:
        symbol_a, symbol_b = ppi.split("*")
        ppi_symbols.append(symbol_a)
        ppi_symbols.append(symbol_b)
    symbol_transcript_map = map_symbols_to_transcripts(ppi_symbols)
    start = time.perf_counter()
    aa_count = await get_sequence_lengths(session, list(symbol_transcript_map.values()))
    finish = time.perf_counter()
    logging.info(f"Made Ensembl API call in {round(finish - start, 2)} second(s)")
    filtered_ppis = []
    for ppi in ppis:
        protein_a, protein_b = ppi.split("*")
        a_transcript = symbol_transcript_map[protein_a]
        b_transcript = symbol_transcript_map[protein_b]
        if aa_count[a_transcript] + aa_count[b_transcript] <= 2000:
            filtered_ppis.append(ppi)
    logging.debug(f"Found a total of {len(filtered_ppis)} ppis from this chunk...")
    return filtered_ppis
