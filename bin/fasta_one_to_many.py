"""
Filter and convert one big FASTA file into 
many smaller files for each transcript.
"""


import argparse
import logging
import os


def process_file(file_path):
    pass


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Input file name')
    parser.add_argument('-o', '--outfolder',
                        type=str, default=None,
                        help='Output folder path (requires absolute path)')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Change logging level from default level to noisiest level')
    parser.add_argument('-l', '--logfile',
                        type=str, default=None,
                        help='Specify the name of the log file')
    args = parser.parse_args()
    return args


def main():
    """Run the command line program."""
    args = parse_command_line()
    outfolder = args.outfolder 
    if args.outfolder is None:
        outfolder = os.path.join(os.getcwd(), args.infile.rsplit('.', 1)[0])
    logfile = args.logfile if args.logfile is not None else 'fasta_one_to_many.log'
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile)
    logging.info(f'Saving files to {outfolder}')
    try:
        process_file(args.infile)
    except FileNotFoundError:
        msg = f'{args.infile} not processed: File does not exist'
        logging.warning(msg)

if __name__ == '__main__':
    main()