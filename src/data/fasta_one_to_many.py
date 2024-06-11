"""
Filter and convert one big FASTA file into 
many smaller files for each transcript.
"""


import argparse
import logging
import os


def _create_files(headers, sequences, outfolder):
    """Write each header sequence pair to its individual FASTA file."""
    for header, sequence in zip(headers, sequences):
        fname = header.split('>')[1].replace('|', '_')
        fname = f'{fname}.fasta'
        file_path = os.path.join(outfolder, fname)
        logging.debug(f'Writing to {file_path}')
        os.makedirs(outfolder, exist_ok=True)
        with open(file_path, 'w') as file:
            file.write(f'{header}\n{sequence}')


def process_file(fname):
    """Split FASTA file into names of smaller 
    FASTA files and their sequences."""
    logging.debug(f'Reading in {fname}...')
    with open(fname, 'r') as reader:
        headers = []
        sequences = []
        current_sequence_pieces = []
        current_header = None
        sequence_available = True
        for line in reader:
            line = line.strip()
            if line.startswith('>'):
                # Process the previous header and sequence
                if current_header:
                    if sequence_available:
                        logging.debug(f'Adding header: {current_header}')
                        headers.append(current_header)
                        current_sequence = '\n'.join(current_sequence_pieces)
                        logging.debug(f'Adding sequence:\n{current_sequence}')
                        sequences.append(current_sequence)
                    current_sequence_pieces = []
                    sequence_available = True
                
                # Start a new header
                current_header = line
            elif line == "Sequence unavailable":
                sequence_available = False
            else:
                if sequence_available:
                    current_sequence_pieces.append(line)
        
        # Process the last header and sequence
        if current_header and sequence_available:
            logging.debug(f'Adding header: {current_header}')
            headers.append(current_header)
            current_sequence = '\n'.join(current_sequence_pieces)
            logging.debug(f'Adding sequence:\n{current_sequence}')
            sequences.append(current_sequence)
    
    return headers, sequences


def parse_command_line(): # pragma: no cover
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


def main(): # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    outfolder = args.outfolder 
    if args.outfolder is None:
        outfolder = os.path.join(os.getcwd(), args.infile.rsplit('.', 1)[0])
    logfile = args.logfile if args.logfile is not None else 'log/fasta_one_to_many.log'
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile)
    try:
        logging.info('Creating header and sequence lists...')
        headers, sequences = process_file(args.infile)
        _create_files(headers, sequences, outfolder)
    except FileNotFoundError:
        msg = f'{args.infile} not processed: File does not exist'
        logging.warning(msg)


if __name__ == '__main__': # pragma: no cover
    main()