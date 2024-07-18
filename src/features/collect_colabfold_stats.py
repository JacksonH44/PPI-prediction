"""
Scan the ColabFold output logs to find the highest pLDDT 
score, the best performing model, and whether there was 
an error in the run.
"""

import argparse
import csv
import logging
import os


def collect_stats(data_dir: str, stats_file: str):
    """
    Scan the data directory and collect stats
    for each new complex created by ColabFold.
    
    Parameters
    ----------
    data_dir : str
        Path to the data directory
    stats_file : str
        Path to the stats file. Assumes stats_file is a CSV
        file with headers symbol,plddt,iptm,best_model,finished
    """
    if os.path.exists(stats_file):
        symbols_needed = os.listdir(data_dir)
        symbols_computed = pd.read_csv(
            stats_file,
            usecols=['symbol']
        )['symbol'].to_list()
        symbols = list(set(symbols_needed) - set(symbols_computed))
    else:
        symbols = os.listdir(data_dir)
        with open(stats_file, 'w') as stats:
            writer = csv.writer(stats)
            writer.writerow(['symbol', 'pLDDT', 'ipTM', 'best_model', 'finished'])
    for symbol in symbols:
        logfile = os.path.join(data_dir, symbol, 'log.txt')
        with open(logfile, 'r') as log:
            lines = log.readlines()
            lines = [line.split(' ', maxsplit=2)[2] for line in lines]
            lines = [line.rstrip('\n') for line in lines]
            finished = False
            plddt = '-1'
            iptm = '-1'
            best_model = 'NA'
            if lines[-1] == 'Done':
                finished = True
                idx = lines.index("reranking models by 'plddt' metric")
                if idx < len(lines) - 1:
                    best = lines[idx + 1]
                    plddt = best.split(' ')[1].split('=')[1]
                    iptm = best.split(' ')[2].split('=')[1]
                    best_model = best.split(' ')[0].split('_')[5]
                    logging.debug(f'pLDDT: {plddt} ipTM: {iptm} best model: {best_model}')
            with open(stats_file, 'a') as stats:
                writer = csv.writer(stats)
                writer.writerow([symbol, plddt, iptm, best_model, finished])


def parse_command_line():  # pragma : no cover
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-d',
        '--data_directory',
        type=str,
        default='data/processed/colabfold',
        help='The path to the directory holding all finished ColabFold outputs',
    )
    parser.add_argument(
        "-s",
        "--stats_file",
        type=str,
        default='data/processed/colabfold_stats.csv',
        help='The path to the file where you wish to write the ColabFold stats',
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

def main():
    """Run the command line program."""
    args = parse_command_line()
    logfile = (
        args.logfile
        if args.logfile is not None
        else os.path.join(os.getcwd(), "logs/collect_colabfold_stats.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    collect_stats(args.data_directory, args.stats_file)


if __name__ == '__main__':
    main()
