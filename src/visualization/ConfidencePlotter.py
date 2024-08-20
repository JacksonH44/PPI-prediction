"""
A class for plotting the difference in AF-M pLDDT
score between positive and negative cases.
"""

import os
import sys

import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from src.visualization.Plotter import Plotter


class ConfidencePlotter(Plotter):
    """
    Plot the difference in AF-M pLDDT scores.

    Usage
    -----
    cp = ConfidencePlotter(
        'data/processed/negative_ppis.csv',
        'data/processed/colabfold_stats.csv'
    )
    cp.plot('reports/figures/confidence.svg')
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
        self._title = "Boxplot of Local Confidence in Predictions"
        self._ylabel = "Local Confidence (AF pLDDT)"

        # Create data object
        stats_df = pd.read_csv(self._stats_file_path, usecols=["symbol", "mean_plddt"])
        neg_conf = []
        pos_conf = []
        for _, row in stats_df.iterrows():
            sym = row["symbol"]
            plddt = row["mean_plddt"]
            if sym in self._negative_ppis:
                neg_conf.append(plddt)
            else:
                pos_conf.append(plddt)

        data = {
            self._ylabel: list(pos_conf) + list(neg_conf),
            "Class": ["Interaction (+)"] * len(pos_conf)
            + ["No Interaction (-)"] * len(neg_conf),
        }
        self._data = pd.DataFrame(data)
        self._pairs = [("Interaction (+)", "No Interaction (-)")]
        self._colour_palette = {
            "Interaction (+)": "#36F1CD",
            "No Interaction (-)": "#D770C7",
        }
