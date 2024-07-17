"""
Utility functions for handling features of
the dataset observations.
"""

import argparse
import csv


def assign_batch(input_file: str, output_file: str, batch_size: int) -> None:
    """
    Assign a batch number to each observation in the input file.

    Parameters
    ----------
    input_file : str
        The path to the file that you want to assign batch numbers to.
        The file should be a CSV with columns symbol,length and the symbol
        should be the name of the MSA that you plan on inputting to ColabFold.
    output_file : str
        The path to the file that the batched input will be saved to.
    batch_size : int
        The number of proteins/complexes per batch.
    """
    new_lines = []
    with open(input_file, "r") as infile:
        reader = csv.reader(infile)
        batch_number = 0
        num_batch_instances = 0
        for row in reader:
            row.append(str(batch_number))
            new_lines.append(row)
            num_batch_instances += 1
            if num_batch_instances == batch_size:
                batch_number += 1
                num_batch_instances = 0

    with open(output_file, "w") as outfile:
        writer = csv.writer(outfile)
        for line in new_lines:
            writer.writerow(line)


def parse_command_line():  # pragma: no cover
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_file",
        type=str,
        help="The path to the file you wish to batch",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        default=None,
        help="Output file path to which you want to save. Not specifying an output file will overwrite the input file.",
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        type=int,
        default=4,
        help="Batch size for running ColabFold",
    )
    args = parser.parse_args()
    return args


def main():  # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    output_file = args.output_file if args.output_file is not None else args.input_file
    assign_batch(args.input_file, output_file, args.batch_size)


if __name__ == "__main__":  # pragma: no cover
    main()
