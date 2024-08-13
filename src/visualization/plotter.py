"""
A class for plotting changes in features of AF-M folded complexes
from the positive and negative set in the dataset.
"""

import pandas as pd

class Plotter:
    """
    An abstract class for plotting the differences in features
    stratified by split (positive or negative).

    Attributes
    ----------
    _negative_ppi_file_path : str
        The path to the file that holds all negative PPIs
    _stats_file_path : str
        The path to the file that holds all stats from the complexes
    """
    def __init__(self, negative_ppi_file_path: str, stats_file_path: str):
        """
        Constructor
        
        Parameters
        ---------
        negative_ppi_file_path : str
            The path to the file that holds all negative PPIs.
        stats_file_path : str
            The path to the file that holds all stats from the complexes.
        """
        self._negative_ppi_file_path = negative_ppi_file_path
        self._stats_file_path = stats_file_path

    def _create_negative_ppis(self):
        """
        Generate a set of negative protein-protein interactions
        used to split the dataset when plotting.
        """
        negative_dataset = pd.read_csv(self._negative_ppi_file_path)
        negative_ppis = set()
        for _, row in negative_dataset.iterrows():
            s = f"{row['gene_symbol_a']}_{row['gene_symbol_b']}"
            negative_ppis.add(s)
        return negative_ppis
