"""
A collection of functions to extract the change
in surface area from a complex.
"""

import argparse
import os
import sys

from contact_map import ContactFrequency  # type: ignore
import mdtraj as md  # type: ignore
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.file_utils import find_pdb_files


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
        # TODO: implement this function
        pass


def parse_command_line():  # pragma: no cover
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "batch_dir",
        type=str,
        help="The path to the batch directory holding files you want to get the surface area for",
    )
    args = parser.parse_args()
    return args


def main():  # pragma: no cover
    """Run the command line program."""
    args = parse_command_line()
    find_delta_surface_areas(args.batch_dir)


if __name__ == "__main__":  # pragma: no cover
    main()
