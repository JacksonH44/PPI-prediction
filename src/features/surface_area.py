"""
A collection of functions to extract the change
in surface area from a complex.
"""

import argparse
import logging
import os
import sys
from typing import Optional

from contact_map import ContactFrequency  # type: ignore
import mdtraj as md  # type: ignore
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_pdb_files


def find_length_split(symbol_name: str, batch_dir: str) -> Optional[tuple[int, int]]:
    """Find the sequence lengths for the two proteins in a protein
    complex."""
    msa_file = f"{batch_dir}/{symbol_name}.msa.a3m"
    try:
        with open(msa_file, "r") as msa:
            firstLine = msa.readline().strip("\n")
            lengths = firstLine.split("\t")[0].lstrip("#").split(",")
            seq_1, seq_2 = int(lengths[0]), int(lengths[1])
            return (seq_1, seq_2)
    except FileNotFoundError:
        logging.error(f"Cannot find processed MSA for {symbol_name} in {batch_dir}")
        return None


def get_contact_map(pdb_file: str) -> pd.DataFrame:
    """Get the contact map for a PDB file. Used script
    from NourHanafi on Github."""
    traj = md.load(pdb_file)
    frame_contacts = ContactFrequency(traj)
    df = frame_contacts.residue_contacts.df

    df = df.fillna(0)
    df = df.where(np.triu(np.ones(df.shape)).astype(bool))
    df = df.stack().reset_index()
    df.columns = ["residue1", "residue2", "weight"]

    df["residue1"] += 1
    df["residue2"] += 1

    df = df.loc[df.residue1.ne(df.residue2)]  # keep non-self rows
    return df


def find_delta_surface_areas(batch_dir: str) -> None:
    """Find the change in surface area for the interaction
    and non-interaction sites for a batch of ColabFold outputs."""
    pdb_files = find_pdb_files(batch_dir)
    for file in pdb_files:
        symbol = file.split("/")[-1].split(".")[0]
        split = find_length_split(symbol, batch_dir)
        if split:
            seq_length_1, seq_length_2 = split[0], split[1]
            logging.debug(f"Sequence lengths for complex {symbol}:")
            logging.debug(
                f"{symbol.split('_')[0]}: {seq_length_1}\n{symbol.split('_')[1]}: {seq_length_2}"
            )


def parse_command_line():  # pragma: no cover
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "batch_dir",
        type=str,
        help="The path to the batch directory holding files you want to get the surface area for",
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
        else os.path.join(os.getcwd(), "logs/surface_area.log")
    )
    os.makedirs(os.path.dirname(logfile), exist_ok=True)
    if not os.path.exists(logfile):
        with open(logfile, "w") as file:
            file.write("")  # Write an empty string to create the file
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=logging_level, filename=logfile, filemode="w")
    find_delta_surface_areas(args.batch_dir)


if __name__ == "__main__":  # pragma: no cover
    main()
