import os
import subprocess
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.features.feature_calculator import FeatureCalculator  # type: ignore


class FrustrationCalculator(FeatureCalculator):
    """
    Calculate the frustration index of a given complex
    """

    def __init__(self, multimer_pdb_path: str, monomer_pdb_path: str):
        """
        Constructor
        """
        super().__init__(multimer_pdb_path, monomer_pdb_path)

    def _get_frst_index(self, pdb_file: str, chain: str) -> list[float]:
        """
        Calculate the average frstIndex across all residues of a PDB file.

        Parameters
        ----------
        pdb_file : str
            The absolute path to the pdb_file you want to compute frustration index for
        chain : str
            The chain you want to compute frustration index for. If the PDB file represents
            a complex, this will be either 'A' or 'B', if the PDB file represents a
            monomer, this will be 'A'

        Returns
        -------
        list[float]
            A list of residue-level frustration index values
        """
        subprocess.run(
            ["Rscript", "src/features/compute_frustration.R", pdb_file, "/tmp", chain],
            stdout=open(os.devnull, "wb"),
        )
        # Get just the file name from the pdb file path
        pdb_file_name = pdb_file.split("/")[-1].rstrip(".pdb")
        output_file = os.path.join(
            "/tmp",
            pdb_file_name + "_" + chain + ".done",
            "FrustrationData",
            pdb_file_name + "_" + chain + ".pdb_singleresidue",
        )
        frst_df = pd.read_csv(output_file, sep=" ", usecols=["FrstIndex"])
        return frst_df["FrstIndex"].to_list()

    def _calculate_residue_metrics(self):
        """
        Calculate the per-residue frustration index for both the
        multimer and monomer
        """
        multimer_chain_no = (
            "A" if self._monomer_symbol == self._complex_symbol.split("_")[0] else "B"
        )
        self._multimer_residue_metrics = self._get_frst_index(
            self._multimer_pdb_path, multimer_chain_no
        )
        self._monomer_residue_metrics = self._get_frst_index(
            self._monomer_pdb_path, "A"
        )


if __name__ == '__main__':
    fc = FrustrationCalculator(
        ('tests/test_data/colabfold/321/SMARCE1_DPF2'
         '.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb'),
        ('tests/test_data/colabfold/monomer/ENST00000528416_DPF2.msa_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb')
    )
    fc.calculate_residue_metrics()
    deltas = fc.calculate_delta_metrics()
    print(deltas)
