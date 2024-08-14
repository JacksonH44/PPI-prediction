"""
A class for plotting the difference in change in
surface area between positive and negative cases.
"""

import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.visualization.Plotter import Plotter


class SAPlotter(Plotter):
    """
    Plot the difference in change in surface areas.
    """

    def __init__(self, negative_ppi_file_path: str, stats_file_path: str):
        """
        Constructor
        """
        super().__init__(negative_ppi_file_path, stats_file_path)

    def _feature_init(self):
        """
        Create features for plotting confidence
        """
        self._title = "Boxplot of Change in Surface Area"
        self._ylabel = "$\Delta$ Surface Area ($\AA$)"

        # Create data object
        stats_df = pd.read_csv(self._stats_file_path, usecols=["symbol", "avg_sa_i_1", "avg_sa_i_2"])
        neg_sa_i_1 = []
        neg_sa_i_2 = []
        pos_sa_i_1 = []
        pos_sa_i_2 = []
        for _, row in stats_df.iterrows():
            sym = row['symbol']
            avg_sa_i_1 = row['avg_sa_i_1']
            avg_sa_i_2 = row['avg_sa_i_2']
            if sym in self._negative_ppis:
                neg_sa_i_1.append(avg_sa_i_1)
                neg_sa_i_2.append(avg_sa_i_2)
            else:
                pos_sa_i_1.append(avg_sa_i_1)
                pos_sa_i_2.append(avg_sa_i_2)
        data = {
            self._ylabel: pos_sa_i_1 + neg_sa_i_1 + pos_sa_i_2 + neg_sa_i_2,
            'Class': ['Cancer Driver (+)'] * len(pos_sa_i_1) + ['Cancer Driver (-)'] * len(neg_sa_i_1) + ['Partner (+)'] * len(pos_sa_i_2) + ['Partner (-)'] * len(neg_sa_i_2)
        }
        self._data = pd.DataFrame(data)
        self._pairs = [('Cancer Driver (+)', 'Cancer Driver (-)'), ('Partner (+)', 'Partner (-)')]
        self._colour_palette = {'Cancer Driver (+)': '#36F1CD', 'Cancer Driver (-)': '#D770C7', 'Partner (+)': '#36F1CD', 'Partner (-)': '#D770C7'}


if __name__ == "__main__":
    sap = SAPlotter(
        "data/processed/negative_ppis.csv", "data/processed/colabfold_stats.csv"
    )
    sap.plot("data/processed/surface_area.svg")
