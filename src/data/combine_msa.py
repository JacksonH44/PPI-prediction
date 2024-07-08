"""
Combine 2 HHblits-generated MSA files into 
1 ColabFold-compliant MSA file.
"""

import argparse
import logging
import os


def extract_header_sequence_pairs(msa_file):
    """Return a hashmap of header-sequence pairs for an
    MSA file."""
    with open(msa_file, 'r') as msa_file:
        lines = msa_file.readlines()
    sequences = {}
    current_header = None
    for line in lines:
        if line.startswith('>'):  # header
            current_header = line.lstrip('>').strip()
        else:  # sequence
            sequences[current_header] = line.strip()
    return sequences


def parse_command_line():  # pragma: no cover
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "first_msa",
        type=str,
        help="The MSA for the first protein in the pair",
    )
    parser.add_argument(
        "second_msa",
        type=str,
        help="The MSA for the second protein in the pair",
    )
    parser.add_argument(
        "output_msa",
        type=str,
        help="The path name to the combined MSA",
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


def main():  # pragma: no cover
    """Run the command-line program."""
    args = parse_command_line()
    logfile = (
        args.logfile
        if args.logfile is not None
        else os.path.join(os.getcwd(), "logs/combine_msa.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    first_msa_sequences = extract_header_sequence_pairs(args.first_msa)
    logging.debug(first_msa_sequences)


if __name__ == '__main__':  # pragma: no cover
    main()