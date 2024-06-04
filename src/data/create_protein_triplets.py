"""
Create dataset observation triplets that have the following
observation structure: (reference protein, alternative isoform, bait protein)
and label structure: [perturbed,conserved]
"""


import argparse
import logging
import os

import pandas as pd


def write_triplet_file(df, outfile):
    """
    Write the triplet dataframe to a CSV.

    Parameters
    ----------
    df : pd.DataFrame
        The observation triplets with labels (perturbation) dataframe.
    outfile : str
        Absolute path to the CSV file that you wish to write the triplets to.
    """
    # Group by T/F for ease of reading the CSV
    df = df.sort_values(by='perturbation', ascending=False)
    df.to_csv(outfile, index=False)


def find_triplets(infile, positive):
    """
    Generate triplets of a reference isoform, an alternative isoform,
    a bait protein, and whether the alternative-bait complex interaction
    changed from the reference-bait complex (the label).

    Parameters
    ----------
    infile : str
        Absolute path to the Excel file holding the protein-protein interactions.
        In the test_data folder, this file is called 'test_ppis.xlsx'.
    positive : Boolean
        A command-line flag specifying whether observations should be restricted 
        to only those where the reference protein interacts with the bait (-p) or
        not. The latter is known as ab initio prediction.
    """
    ppi_df = pd.read_excel(infile, sheet_name='2B-Isoform PPIs', 
                           usecols=['Gene_Symbol','Isoform_ID', 'Category', 
                                    'Interactor_ID', 'Interaction_Found'])
    if positive:
        reference_isoforms = ppi_df[(ppi_df['Category'] == 'reference') & (ppi_df['Interaction_Found'] == 'positive')]
    else:
        reference_isoforms = ppi_df[(ppi_df['Category'] == 'reference')]
    logging.debug(f'Reference proteins:\n{reference_isoforms}')
    observations = []

    for gene_symbol in reference_isoforms['Gene_Symbol'].unique():
        isoform_id = f'{gene_symbol}_1'
        logging.debug(f'Isoform: {isoform_id}, Gene symbol: {gene_symbol}')
        # Remove alternative isoforms that don't interact with the bait protein
        alternative_isoforms = ppi_df[(ppi_df['Gene_Symbol'] == gene_symbol) & (ppi_df['Category'] == 'alternative') 
                                      & (pd.notna(ppi_df['Interaction_Found']))]
        logging.debug(f'Alternative isoforms for {isoform_id}:\n{alternative_isoforms}')

        for _, alt_row in alternative_isoforms.iterrows():
            ref_isoform = reference_isoforms[(reference_isoforms['Gene_Symbol'] == gene_symbol) 
                                             & (reference_isoforms['Interactor_ID'] == alt_row['Interactor_ID'])]
            if ref_isoform.empty:
                continue
            logging.debug(f'Alternative isoform:\n{alt_row}\nReference isoform:\n{ref_isoform}')
            # Create 4-tuples for the isoform-bait interaction
            ref_id = ref_isoform['Isoform_ID'].values[0]
            alt_id = alt_row['Isoform_ID']
            bait_id = alt_row['Interactor_ID']
            perturbation = ref_isoform['Interaction_Found'].values[0] != alt_row['Interaction_Found']
            logging.debug(f'Adding observation ({ref_id}, {alt_id}, {bait_id}, {perturbation})')
            observations.append([ref_id, alt_id, bait_id, perturbation])
    observation_df = pd.DataFrame(observations, columns=['ref_ID', 'alt_ID', 'bait_ID', 'perturbation'])
    logging.debug(f'Observations:\n{observation_df}')
    return observation_df


def parse_command_line():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', type=str,
                        help='Absolute path to the protein-protein interaction Excel file.')
    parser.add_argument('outfile', type=str,
                        help='Absolute path of the CSV file to write to')
    parser.add_argument('-p', '--positive',
                        action='store_true',
                        help='''Include this flag if you want the dataset to include only 
                        reference isoforms that interact with the bait protein''')
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
    logfile = args.logfile if args.logfile is not None else os.path.join(os.getcwd(), "logs/create_protein_triplets.log")
    if not os.path.exists(logfile):
        with open(logfile, 'w') as file:
            file.write("") # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, filename=logfile, filemode='w')
    logging.debug('Creating triplets...')
    observations = find_triplets(args.infile, args.positive)
    logging.debug(f'Writing to triplet file in location {args.outfile}')
    write_triplet_file(observations, args.outfile)


if __name__ == '__main__':
    main()