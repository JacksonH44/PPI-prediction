"""
Scan the ColabFold output logs to find the highest pLDDT
score, the best performing model, and whether there was
an error in the run.
"""

import argparse
import csv
import json
import logging
import os
import sys
import time
import re

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join("..", "..")))
from src.features.file_utils import find_all_complexes
from src.features.find_stats import find_stats


def get_colabfold_metrics(symbol: str, data_dir: str) -> list[str]:
    """
    Return the max predicted aligned error (PAE), the pLDDT score, ipTM score,
    and pTM score for a given complex.

    Parameters
    ----------
    symbol : str
        The complex symbol
    data_dir : str
        The path to the directory containing the scores for the complex

    Returns
    -------
    str_features : list[str]
        A 10 element list representing 10 features: max PAE, average, max, and
        min for each of pLDDT, pTM, and ipTM.
    """
    ptm_regex = rf"^{symbol}.msa_scores_rank_00\d+_alphafold2_multimer_v3_model_\d+_seed_000.json"
    ptm_scores = [file for file in os.listdir(data_dir) if re.search(ptm_regex, file)]
    max_ptm, max_iptm, max_plddt = float("-inf"), float("-inf"), float("-inf")
    min_ptm, min_iptm, min_plddt = float("inf"), float("inf"), float("inf")
    avg_ptm, avg_iptm, avg_plddt = 0.0, 0.0, 0.0
    for ptm_score_file in ptm_scores:
        with open(os.path.join(data_dir, ptm_score_file), "r") as ptm_json:
            score = json.load(ptm_json)
            ptm = score["ptm"]
            iptm = score["iptm"]
            raw_plddt = score["plddt"]
            plddt = sum(raw_plddt) / len(raw_plddt)
            max_ptm = max(max_ptm, ptm)
            max_iptm = max(max_iptm, iptm)
            max_plddt = max(max_plddt, plddt)
            min_ptm = min(min_ptm, ptm)
            min_iptm = min(min_iptm, iptm)
            min_plddt = min(min_plddt, plddt)
            avg_ptm += ptm
            avg_iptm += iptm
            avg_plddt += plddt
    avg_ptm = round(avg_ptm / len(ptm_scores) if avg_ptm != 0.0 else 0.0, 2)
    avg_iptm = round(avg_iptm / len(ptm_scores) if avg_iptm != 0.0 else 0.0, 2)
    avg_plddt = round(avg_plddt / len(ptm_scores) if avg_plddt != 0.0 else 0.0, 2)
    max_plddt = round(max_plddt, 4)
    min_plddt = round(min_plddt, 4)

    # Get predicted aligned error score
    pae_regex = rf"^{symbol}.msa_predicted_aligned_error_v1.json"
    pae_file = [file for file in os.listdir(data_dir) if re.search(pae_regex, file)]
    assert len(pae_file) == 1, f"PAE file for {symbol} is missing"
    with open(os.path.join(data_dir, pae_file[0]), "r") as pae_json:
        score = json.load(pae_json)
        pae = score["max_predicted_aligned_error"]

    float_features = [
        pae,
        avg_plddt,
        max_plddt,
        min_plddt,
        avg_ptm,
        max_ptm,
        min_ptm,
        avg_iptm,
        max_iptm,
        min_iptm,
    ]
    str_features = [str(feat) for feat in float_features]
    return str_features


def collect_stats(data_dir: str, monomer_dir: str, stats_file: str):
    """
    Scan the data directory and collect stats
    for each new complex created by ColabFold.

    Parameters
    ----------
    data_dir : str
        Path to the data directory
    monomer_dir : str
        Path to the monomer file directory
    stats_file : str
        Path to the stats file. Assumes stats_file is a CSV
        file with headers symbol,plddt,iptm,best_model
    """
    symbols = find_all_complexes(data_dir)
    logging.debug(f"Found complexes:\n{symbols}")
    with open(stats_file, "w") as stats:
        writer = csv.writer(stats)
        header = [
            "symbol",
            "max_pae",
            "mean_plddt",
            "max_plddt",
            "min_plddt",
            "mean_ptm",
            "max_ptm",
            "min_ptm",
            "mean_iptm",
            "max_iptm",
            "min_iptm",
            "min_sa_i_1",
            "min_sa_i_2",
            "max_sa_i_1",
            "max_sa_i_2",
            "avg_sa_i_1",
            "avg_sa_i_2",
            "min_sa_ni_1",
            "min_sa_ni_2",
            "max_sa_ni_1",
            "max_sa_ni_2",
            "avg_sa_ni_1",
            "avg_sa_ni_2",
            "min_fi_i_1",
            "min_fi_i_2",
            "max_fi_i_1",
            "max_fi_i_2",
            "avg_fi_i_1",
            "avg_fi_i_2",
            "min_fi_ni_1",
            "min_fi_ni_2",
            "max_fi_ni_1",
            "max_fi_ni_2",
            "avg_fi_ni_1",
            "avg_fi_ni_2",
        ]
        logging.debug("Writing header...")
        writer.writerow(header)
        for symbol in symbols:
            try:
                features = [symbol]
                colabfold_metrics = get_colabfold_metrics(symbol, data_dir)
                features += colabfold_metrics
                logging.debug(f"Finding surface area for {symbol}")
                features += find_stats(symbol, data_dir, monomer_dir, "surface_area", 5)
                features += find_stats(symbol, data_dir, monomer_dir, "frustration", 5)
                writer.writerow(features)
            except AssertionError as e:
                print(f"Could not find features for {symbol}\n{e}")


def parse_command_line():  # pragma : no cover
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-d",
        "--data_directory",
        type=str,
        default=os.path.join("tests", "test_data", "colabfold", "0"),
        help="The path to the directory holding all finished ColabFold outputs",
    )
    parser.add_argument(
        "-m",
        "--monomer_directory",
        type=str,
        default=os.path.join("tests", "test_data", "colabfold", "monomer"),
        help="The path to the directory holding all monomer ColabFold outputs",
    )
    parser.add_argument(
        "-s",
        "--stats_file",
        type=str,
        default=os.path.join("data", "processed", "colabfold_stats.csv"),
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
        else os.path.join(
            os.getcwd(), os.path.join("logs", "collect_colabfold_stats.log")
        )
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    start = time.perf_counter()  # Time the function call for debugging
    collect_stats(args.data_directory, args.monomer_directory, args.stats_file)
    finish = time.perf_counter()
    logging.info(f"Finished curation in {round(finish - start, 2)} second(s)")


if __name__ == "__main__":  # pragma: no cover
    main()
