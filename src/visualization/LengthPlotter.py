"""
A class for plotting the difference in lengths
between a positive and negative case.
"""

import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join("..", "..")))
from src.visualization.Plotter import Plotter


class LengthPlotter(Plotter):
    """
    Plot the difference in lengths

    Attributes
    ----------
    _length_file_path : str
        The name of the path to which the complex lengths are stored.

    Usage
    -----
    lp = LengthPlotter(
        os.path.join('data', 'processed', 'negative_ppis.csv'),
        os.path.join('data', 'processed', 'colabfold_stats.csv')
    )
    lp.plot(os.path.join('reports', 'figures', 'length.svg'))
    """

    def __init__(
        self, negative_ppi_file_path: str, stats_file_path: str, length_file_path: str
    ):
        """
        Constructor

        Parameters
        ----------
        length_file_path : str
            The name of the path to which the complex lengths are stored
        """
        self._length_file_path = length_file_path
        super().__init__(negative_ppi_file_path, stats_file_path)

    def _feature_init(self):
        """
        Create features for plotting lengths
        """
        self._title = "Boxplot of Amino Acid Lengths"
        self._ylabel = "Amino Acid Length"
        lengths_df = pd.read_csv(self._length_file_path, usecols=["symbol", "length"])

        # Create data object
        stats_df = pd.read_csv(self._stats_file_path, usecols=["symbol"])
        neg_lengths = []
        pos_lengths = []
        for _, row in stats_df.iterrows():
            sym = row["symbol"]
            length = lengths_df[lengths_df["symbol"] == row["symbol"]]["length"].values[
                0
            ]
            if sym in self._negative_ppis:
                neg_lengths.append(length)
            else:
                pos_lengths.append(length)

        data = {
            self._ylabel: list(pos_lengths) + list(neg_lengths),
            "Class": ["Interaction (+)"] * len(pos_lengths)
            + ["No Interaction (-)"] * len(neg_lengths),
        }
        self._data = pd.DataFrame(data)
        self._pairs = [("Interaction (+)", "No Interaction (-)")]
        self._colour_palette = {
            "Interaction (+)": "#36F1CD",
            "No Interaction (-)": "#D770C7",
        }
