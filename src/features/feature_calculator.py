from abc import ABC, abstractmethod
import os
from typing import Optional

import numpy as np
from contact_map import ContactFrequency  # type: ignore
import mdtraj as md  # type: ignore


class FeatureCalculator(ABC):
    """
    A class that represents a feature calculation

    Attributes
    ----------
    _multimer_pdb_path : str
        The absolute path to the multimer PDB file
    _monomer_pdb_path : str
        The absolute path to the monomer PDB file
    _complex_symbol : str
        The EntrezGene symbol of the two genes in the complex
        joined with a '_'
    _monomer_symbol : str
        The EntrezGene symbol of the gene in the monomer file
    _contact_map : str
        A dataframe representing a contact map for the complex,
        if a row contains residue X and residue Y, then X (from chain A)
        is in contact with Y (from chain B)
    _mask_map : dict[str, list[bool]]
        A map of complex, mask pairs, where the mask[i] is True if
        residue i + 1 is in the interaction site of the complex, and
        mask[i] is False if residue i + 1 is not part of the interaction
        site
    _residue_metrics : list[float]
        The metrics for each residue in the complex
    """

    def __init__(self, multimer_pdb_path: str, monomer_pdb_path: str):
        """
        Constructor

        Parameters
        ----------
        multimer_pdb_path : str
            The absolute path to the multimer PDB file
        monomer_pdb_path : str
            The absolute path to the monomer PDB file
        """
        self._multimer_pdb_path = multimer_pdb_path
        self._monomer_pdb_path = monomer_pdb_path
        self._complex_symbol = os.path.basename(self._multimer_pdb_path).split(".")[0]
        self._monomer_symbol = (
            os.path.basename(self._monomer_pdb_path).split(".")[0].split("_")[-1]
        )
        self._create_length_split()
        self._create_contact_map()
        self._create_interaction_site()
        self._multimer_residue_metrics = None
        self._monomer_residue_metrics = None
        self._multimer_interaction_site = None
        self._multimer_non_interaction_site = None
        self._monomer_interaction_site = None
        self._monomer_non_interaction_site = None

    @property
    def complex_symbol(self):
        return self._complex_symbol

    @property
    def monomer_symbol(self):
        return self._monomer_symbol

    def _create_length_split(self):
        """
        Find the sequence lengths for the two proteins in a protein
        complex. NOTE: Assumes that the MSA file for the complex is stored in
        the same directory the PDB output file is, which is the case for
        ColabFold outputs.

        Returns
        -------
        tuple[int, int], optional
            A tuple of the cancer driver sequence length, and the interactor
            sequence length
        """
        msa_file = os.path.join(
            os.path.dirname(self._multimer_pdb_path), f'{self._complex_symbol}.msa.a3m'
        )
        try:
            with open(msa_file, "r") as msa:
                firstLine = msa.readline().strip("\n")
                lengths = firstLine.split("\t")[0].lstrip("#").split(",")
                self._seq_1_length, self._seq_2_length = int(lengths[0]), int(
                    lengths[1]
                )
        except FileNotFoundError:
            print(f"Cannot find processed MSA for {self._complex_symbol}")
            self._seq_1_length, self._seq_2_length = None, None

    def _create_contact_map(self):
        """
        Get the contact map for a PDB file. Used script
        from NourHanafi on Github.

        Parameters
        ----------
        pdb_file : str
            The path to the PDB file for a complex to obtain the
            contact map for

        Returns
        -------
        pd.DataFrame
            A dataframe containing the residue numbers that are
            in contact
        """
        traj = md.load(self._multimer_pdb_path)
        frame_contacts = ContactFrequency(traj)
        df = frame_contacts.residue_contacts.df

        df = df.fillna(0)
        df = df.where(np.triu(np.ones(df.shape)).astype(bool))
        df = df.stack().reset_index()
        df.columns = ["residue1", "residue2", "weight"]

        df["residue1"] += 1
        df["residue2"] += 1

        df = df.loc[df.residue1.ne(df.residue2)]  # keep non-self rows
        df = df[df["weight"] == 1.0]
        self._contact_map = df[["residue1", "residue2"]]

    def _create_interaction_site(self):
        """
        Find the interaction sites for the two proteins present in the complex. It considers
        cases where one residue in the contact map is on one protein, while the other residue
        is on the other protein. Indices returned in the mask are relative to the single protein
        (e.g. if an interacting residue is residue #214 in the complex but #89 in the single
        protein, the number 89 will be returned).
        """
        protein_1, protein_2 = (
            self._complex_symbol.split("_")[0],
            self._complex_symbol.split("_")[1],
        )
        mask_1 = [False] * self._seq_1_length
        mask_2 = [False] * self._seq_2_length
        for row in self._contact_map.itertuples(index=False):
            residue_1, residue_2 = row[0], row[1]
            # NOTE: Indices of the masks are one less than the residues because
            # residues are 1-indexed and masks are 0-indexed
            # residue 1 is on protein 1 and residue 2 on protein 2
            if residue_1 <= self._seq_1_length and residue_2 > self._seq_1_length:
                mask_1[residue_1 - 1] = True
                mask_2[residue_2 - self._seq_1_length - 1] = True
            # residue 1 is on protein 2 and residue 2 on protein 1
            if residue_1 > self._seq_1_length and residue_2 <= self._seq_1_length:
                mask_2[residue_1 - self._seq_1_length - 1] = True
                mask_1[residue_2 - 1] = True
        self._mask_map = {}
        self._mask_map[protein_1] = mask_1
        self._mask_map[protein_2] = mask_2

    def _apply_residue_mask(
        self, residue_metrics: list[float]
    ) -> Optional[tuple[list[float], list[float]]]:
        """
        Get a list of residues and a mask as input, and return two lists of residues -
        one for the interaction site, and one for the non-interaction site
        """
        if (
            self._multimer_residue_metrics is None
            or self._monomer_residue_metrics is None
        ):
            print(
                "Calculate the residue-level metric before applying the interaction site mask"
            )
            return None
        symbol_mask = self._mask_map[self._monomer_symbol]
        interaction_site = [
            residue_metrics[i] for i in range(len(residue_metrics)) if symbol_mask[i]
        ]
        non_interaction_site = [
            residue_metrics[i]
            for i in range(len(residue_metrics))
            if not symbol_mask[i]
        ]
        return interaction_site, non_interaction_site

    @abstractmethod
    def _calculate_residue_metrics(self):
        """
        This should calculate residue level metrics, and assign values
        to self._multimer_residue_metrics and self._monomer_residue_metrics
        """

    def calculate_residue_metrics(self):
        """
        Calculate per-residue level features, then split the features
        into interaction site and non-interaction site for both the multimer
        and monomer
        """
        self._calculate_residue_metrics()
        self._multimer_interaction_site, self._multimer_non_interaction_site = (
            self._apply_residue_mask(self._multimer_residue_metrics)
        )
        self._monomer_interaction_site, self._monomer_non_interaction_site = (
            self._apply_residue_mask(self._monomer_residue_metrics)
        )

    def calculate_delta_metrics(self) -> tuple[float, float]:
        """
        Return the average difference in the multimer metric and the monomer metric

        Returns
        -------
        tuple[float, float]
            The difference in average metric between multimer and monomer for the interaction site,
            non-interaction site
        """
        assert (
            self._monomer_interaction_site is not None
            and self._multimer_interaction_site is not None
            and self._multimer_non_interaction_site is not None
            and self._monomer_non_interaction_site is not None
        ), "Interaction site splits should not have None type"
        assert len(self._monomer_interaction_site) == len(
            self._multimer_interaction_site
        ), "The lengths of monomer and multimer residue interaction site should be the same"
        assert len(self._monomer_non_interaction_site) == len(
            self._multimer_non_interaction_site
        ), "The lengths of monomer and multimer residue interaction site should be the same"
        interaction_delta = round(
            (sum(self._multimer_interaction_site) - sum(self._monomer_interaction_site))
            / len(self._monomer_interaction_site) if len(self._monomer_interaction_site) != 0 else 0,
            4,
        )
        non_interaction_delta = round(
            (
                sum(self._multimer_non_interaction_site)
                - sum(self._monomer_non_interaction_site)
            )
            / len(self._monomer_non_interaction_site),
            4,
        )
        return (interaction_delta, non_interaction_delta)
