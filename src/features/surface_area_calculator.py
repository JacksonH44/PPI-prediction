import os
import sys
from Bio.PDB import PDBParser, SASA  # type: ignore

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join("..", "..")))
from src.features.feature_calculator import FeatureCalculator  # type: ignore


class SurfaceAreaCalculator(FeatureCalculator):
    """
    Calculate the surface area for a given complex

    Usage
    -----
    sac = SurfaceAreaCalculator(
        os.path.join(
            'tests',
            'test_data',
            'colabfold',
            '321',
            'SMARCE1_DPF2.msa_unrelaxed_rank_001_alphafold2_multimer_v3_model_4_seed_000.pdb'
        ),
        os.path.join(
            'tests',
            'test_data',
            'colabfold',
            'monomer',
            'ENST00000528416_DPF2.msa_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb'
        )
    )
    sac.calculate_residue_metrics()
    deltas = sac.calculate_delta_metrics()
    """

    def __init__(self, multimer_pdb_path: str, monomer_pdb_path: str):
        """
        Constructor
        """
        super().__init__(multimer_pdb_path, monomer_pdb_path)

    def _calculate_surface_area_struct(self, pdb_path: str):
        """
        Calculates surface areas for each residue of a structure in a PDB file

        Parameters
        ----------
        pdb_path : str
            The path to the pdb file to create a surface area structure object from

        Returns
        -------
        sr : Structure
            A structure object representing the residue-level computed surface areas
            for the PDB file
        """
        p = PDBParser(QUIET=1)
        symbol = os.path.basename(pdb_path).split(".")[0]
        struct = p.get_structure(symbol, pdb_path)
        sr = SASA.ShrakeRupley()
        sr.compute(struct, level="R")
        return struct

    def _extract_residues(self, chain) -> list[float]:
        """
        Return a list with the same length as the input sequence of surface areas,
        one entry for each sequence. The model index and chain should already be specified
        (e.g. usage: extract_residues(struct[0]['A']))

        Parameters
        ----------
        chain : Chain
            The input chain for which to extract residue-level surface area calculations from

        Returns
        -------
        surface_areas : list[float]
            A list of surface areas the same length as the number of residues in the chain rounded
            to 4 significant digits
        """
        surface_areas = [round(chain[i].sasa, 4) for i in range(1, len(chain) + 1)]
        return surface_areas

    def _calculate_residue_metrics(self):
        """
        Calculate the per-residue surface area for both the
        multimer and monomer
        """
        multimer_chain_no = (
            "A" if self._monomer_symbol == self._complex_symbol.split("_")[0] else "B"
        )
        multimer_struct = self._calculate_surface_area_struct(self._multimer_pdb_path)
        self._multimer_residue_metrics = self._extract_residues(
            multimer_struct[0][multimer_chain_no]
        )
        monomer_struct = self._calculate_surface_area_struct(self._monomer_pdb_path)
        self._monomer_residue_metrics = self._extract_residues(monomer_struct[0]["A"])
