"""
Combine 2 HHblits-generated MSA files into
1 ColabFold-compliant MSA file.
"""

import argparse
import logging
import os


def write_combined_a3m(
    first_msa_sequences: dict[str, str],
    second_msa_sequences: dict[str, str],
    output_msa: str,
) -> None:
    """Write the header-sequence pairs as one MSA file in
    accordance with the example shown in ColabFold issue
    #76."""
    first_length = len(list(first_msa_sequences.values())[0])
    second_length = len(list(second_msa_sequences.values())[0])
    if os.path.isabs(output_msa):
        os.makedirs(os.path.dirname(output_msa), exist_ok=True)
    else:
        os.makedirs(
            os.path.dirname(os.path.join(os.getcwd(), output_msa)), exist_ok=True
        )
    with open(output_msa, "w") as combined_msa:
        # Write MSA header
        combined_msa.write(f"#{first_length},{second_length}\t1,1\n")
        # Write paired sequence portion of MSA
        for (first_header, first_sequence), (second_header, second_sequence) in zip(
            first_msa_sequences.items(), second_msa_sequences.items()
        ):
            combined_msa.write(
                f">{first_header}\t{second_header}\n{first_sequence}{second_sequence}\n"
            )
        # Write unpaired portion for first MSA
        for header, sequence in first_msa_sequences.items():
            combined_msa.write(f'>{header}\n{sequence}{"-" * second_length}\n')
        # Write unpaired portion for second MSA
        for header, sequence in second_msa_sequences.items():
            combined_msa.write(f'>{header}\n{"-" * first_length}{sequence}\n')


def extract_header_sequence_pairs(msa_file) -> dict[str, str]:
    """Return a hashmap of header-sequence pairs for an
    MSA file."""
    with open(msa_file, "r") as msa_file:
        lines = msa_file.readlines()
    sequences: dict[str, str] = {}
    current_header = ""
    for line in lines:
        if line.startswith(">"):  # header
            current_header = line.lstrip(">").split(" ")[0].strip()
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
    second_msa_sequences = extract_header_sequence_pairs(args.second_msa)
    write_combined_a3m(first_msa_sequences, second_msa_sequences, args.output_msa)


if __name__ == "__main__":  # pragma: no cover
    main()
