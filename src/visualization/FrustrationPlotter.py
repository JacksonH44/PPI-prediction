"""
A class for plotting the difference in change in
frustration index between positive and negative cases.
"""

import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.visualization.Plotter import Plotter


class FrustrationPlotter(Plotter):
    """
    Plot the difference in change in frustration indices.
    """

    def __init__(self, negative_ppi_file_path: str, stats_file_path: str):
        """
        Constructor
        """
        super().__init__(negative_ppi_file_path, stats_file_path)

    def _feature_init(self):
        """
        Create features for plotting frustration index
        """
        self._title = "Boxplot of Change in Energetic Frustration"
        self._ylabel = "$\Delta$ Frustration Index"

        # Create data object
        stats_df = pd.read_csv(self._stats_file_path, usecols=["symbol", 'avg_fi_i_1', 'avg_fi_i_2'])
        neg_fi_i_1 = []
        neg_fi_i_2 = []
        pos_fi_i_1 = []
        pos_fi_i_2 = []
        for _, row in stats_df.iterrows():
            sym = row['symbol']
            avg_fi_i_1 = row['avg_fi_i_1']
            avg_fi_i_2 = row['avg_fi_i_2']
            if sym in self._negative_ppis:
                neg_fi_i_1.append(avg_fi_i_1)
                neg_fi_i_2.append(avg_fi_i_2)
            else:
                pos_fi_i_1.append(avg_fi_i_1)
                pos_fi_i_2.append(avg_fi_i_2)
        data = {
            self._ylabel: pos_fi_i_1 + neg_fi_i_1 + pos_fi_i_2 + neg_fi_i_2,
            'Class': ['Cancer Driver (+)'] * len(pos_fi_i_1) + ['Cancer Driver (-)'] * len(neg_fi_i_1) + ['Partner (+)'] * len(pos_fi_i_2) + ['Partner (-)'] * len(neg_fi_i_2)
        }
        self._data = pd.DataFrame(data)
        self._pairs = [('Cancer Driver (+)', 'Cancer Driver (-)'), ('Partner (+)', 'Partner (-)')]
        self._colour_palette = {'Cancer Driver (+)': '#36F1CD', 'Cancer Driver (-)': '#D770C7', 'Partner (+)': '#36F1CD', 'Partner (-)': '#D770C7'}


if __name__ == "__main__":
    sap = FrustrationPlotter(
        "data/processed/negative_ppis.csv", "data/processed/colabfold_stats.csv"
    )
    sap.plot("data/processed/frustration.svg")
