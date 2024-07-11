"""
Get a report of all sequence/complex lengths, and write
them to a file. Additionally, generate a figure of the 
distribution of protein sequence lengths.
"""

import asyncio
import argparse
import logging
import math
import os
import sys

import aiohttp
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.data.bio_apis import get_sequence_lengths
from src.data.data_processing import (
    chunk_input_genes,
    find_unique_genes,
    map_symbols_to_transcripts
)


def plot_lengths_distribution(lengths: list[str], save_path: str) -> None:
    """Plot a histogram of protein lengths."""
    num_bins = int(math.sqrt(len(lengths)))
    plt.hist(lengths, bins=num_bins)
    plt.xlabel('Sequence lengths')
    plt.ylabel('Frequency')
    plt.title('Distribution of Protein Sequence Lengths')
    plt.savefig(save_path)


def parse_command_line():  # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-p',
        '--positive',
        type=str,
        default='data/processed/positive_ppis.csv',
        help='Path to the positive PPIs dataset'
    )
    parser.add_argument(
        '-n',
        '--negative',
        type=str,
        default='data/processed/negative_ppis.csv',
        help='Path to the negative PPIs dataset'
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        default="data/interim/sequence_lengths.csv",
        help="Path to output file to store sequence lengths",
    )
    parser.add_argument(
        "-i",
        "--visualization",
        type=str,
        default="data/processed/length_distribution.png",
        help="Path to output file of lengths histogram",
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
        else os.path.join(os.getcwd(), "logs/find_sequence_lengths.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    all_genes = find_unique_genes([args.positive, args.negative])
    all_transcripts = map_symbols_to_transcripts(list(all_genes))
    chunked_transcripts = chunk_input_genes(list(all_transcripts.values()), chunk_size=200)
    async with aiohttp.ClientSession() as session:
        tasks = [get_sequence_lengths(session, chunk) for chunk in chunked_transcripts]
        mapped_sequences = await asyncio.gather(*tasks, return_exceptions=True)
    if isinstance(mapped_sequences, BaseException):
        logging.warning('Error in async API call. Exiting...')
        return
    transcript_map = {
        transcript: length
        for hashmap in mapped_sequences
        for transcript, length in hashmap.items()  # type: ignore
    }
    sorted_transcript_map = {transcript: length for transcript, length in sorted(transcript_map.items(), key = lambda item: item[1])}
    logging.debug('Plotting histogram of sequence lengths')
    plot_lengths_distribution(list(sorted_transcript_map.values()), args.visualization)


if __name__ == '__main__':  # pragma: no cover
    asyncio.run(main())
