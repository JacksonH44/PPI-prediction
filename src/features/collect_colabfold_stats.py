"""
Scan the ColabFold output logs to find the highest pLDDT
score, the best performing model, and whether there was
an error in the run.
"""

import argparse
import csv
import logging
import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_all_complexes


def _extract_metrics_from_line(line: str) -> tuple[float, float, float]:
    """
    Take a line of the form:

    rank_001_alphafold2_multimer_v3_model_4_seed_000 pLDDT=70.2 pTM=0.663 ipTM=0.135

    and return pLDDT, pTM, and ipTM
    """
    plddt = line.split(" ")[1].split("=")[1]
    ptm = line.split(" ")[2].split("=")[1]
    iptm = line.split(" ")[3].split("=")[1]
    return float(plddt), float(ptm), float(iptm)


def get_colabfold_metrics(symbol: str, lines: list[str]) -> list[str]:
    """
    Return the best model, pLDDT score, and ipTM score from the ColabFold log
    output file.

    Parameters
    ----------
    symbol : str
        The complex symbol
    lines : list[str]
        The processed lines of the logfile ready for metric extraction

    Returns
    -------
    list[str]
        A 13 element list representing the best model (first element in list)
        and 12 features: best model, average, max, and min for each of pLDDT,
        pTM, and ipTM.
    """
    # Find the start and the end of the log file for this symbol
    # Line should say something like 'Query 1/4 CDKN2C_CD24.msa (length 248)'
    query_start = [line for line in lines if symbol in line][0]
    query_start_index = lines.index(query_start)
    query_end = [
        lines[i]
        for i in range(query_start_index + 1, len(lines))
        if "Query" in lines[i] or "Done" in lines[i]
    ][0]
    query_end_index = lines.index(query_end)
    lines = lines[query_start_index:query_end_index]
    stats_lines = lines[-5:]
    best_model = stats_lines[0].split(" ")[0].split("_")[6]
    best_model_plddt, best_model_ptm, best_model_iptm = _extract_metrics_from_line(
        stats_lines[0]
    )
    plddt_total = 0.0
    ptm_total = 0.0
    iptm_total = 0.0
    max_plddt = 0.0
    min_plddt = 100.0
    max_ptm = 0.0
    min_ptm = 1.0
    max_iptm = 0.0
    min_iptm = 1.0
    for line in stats_lines:
        plddt, ptm, iptm = _extract_metrics_from_line(line)
        plddt_total += plddt
        max_plddt = max(max_plddt, plddt)
        min_plddt = min(min_plddt, plddt)
        ptm_total += ptm
        max_ptm = max(max_ptm, ptm)
        min_ptm = min(min_ptm, ptm)
        iptm_total += iptm
        max_iptm = max(max_iptm, iptm)
        min_iptm = min(min_iptm, iptm)
    mean_plddt = round(plddt_total / 5, 4)
    mean_ptm = round(ptm_total / 5, 4)
    mean_iptm = round(iptm_total / 5, 4)

    features = [
        best_model,
        best_model_plddt,
        mean_plddt,
        max_plddt,
        min_plddt,
        best_model_ptm,
        mean_ptm,
        max_ptm,
        min_ptm,
        best_model_iptm,
        mean_iptm,
        max_iptm,
        min_iptm,
    ]
    features_str = [str(feat) for feat in features]
    return features_str


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
        file with headers symbol,plddt,iptm,best_model
    """
    if os.path.exists(stats_file):
        symbols_needed = find_all_complexes(data_dir)
        symbols_computed = pd.read_csv(stats_file, usecols=["symbol"])[
            "symbol"
        ].to_list()
        symbols = list(set(symbols_needed) - set(symbols_computed))
    else:
        symbols = find_all_complexes(data_dir)
        with open(stats_file, "w") as stats:
            writer = csv.writer(stats)
            header = [
                "symbol",
                "best_model",
                "best_model_plddt",
                "mean_plddt",
                "max_plddt",
                "min_plddt",
                "best_model_ptm",
                "mean_ptm",
                "max_ptm",
                "min_ptm",
                "best_model_iptm",
                "mean_iptm",
                "max_iptm",
                "min_iptm",
            ]
            writer.writerow(header)
    logfile = os.path.join(data_dir, "log.txt")
    with open(logfile, "r") as log:
        lines = log.readlines()
        lines = [line.split(" ", maxsplit=2)[2] for line in lines]
        lines = [line.rstrip("\n") for line in lines]
        for symbol in symbols:
            features = []
            features.append(symbol)
            features += get_colabfold_metrics(symbol, lines)
            with open(stats_file, "a") as stats:
                writer = csv.writer(stats)
                writer.writerow(features)


def parse_command_line():  # pragma : no cover
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-d",
        "--data_directory",
        type=str,
        default="data/processed/colabfold/0",
        help="The path to the directory holding all finished ColabFold outputs",
    )
    parser.add_argument(
        "-s",
        "--stats_file",
        type=str,
        default="data/processed/colabfold_stats.csv",
        help="The path to the file where you wish to write the ColabFold stats",
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


if __name__ == "__main__":  # pragma: no cover
    main()
