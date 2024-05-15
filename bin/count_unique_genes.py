"""
Count the number of unique genes in a TSV
file of interacting protein-coding genes.
"""


import argparse
import csv


def process_file(file_path):
    """Count the number of unique genes in the file."""
    unique_genes = set()

    with open(file_path, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            gene1, gene2 = row
            unique_genes.add(gene1)
            unique_genes.add(gene2)

    return len(unique_genes)


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Input file name')
    args = parser.parse_args()
    return args


def main():
    """Run the command line program."""
    args = parse_command_line()
    num_unique_genes = process_file(args.infile)
    print(f'There are {num_unique_genes} unique genes in the file')


if __name__ == '__main__':
    main()