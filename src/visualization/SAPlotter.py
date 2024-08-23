"""
A class for plotting the difference in change in
surface area between positive and negative cases.
"""

from collections import defaultdict
import os
import sys

import pandas as pd
import matplotlib.pyplot as plt  # type: ignore

sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join("..", "..")))
from src.visualization.Plotter import Plotter


class SAPlotter(Plotter):
    """
    Plot the difference in change in surface areas.

    Usage
    -----
    sap = SAPlotter(
        os.path.join('data', 'processed', 'negative_ppis.csv'),
        os.path.join('data', 'processed', 'colabfold_stats.csv')
    )
    sap.plot(os.path.join('reports', 'figures', 'confidence.svg'))
    sap.plot_biggest_sa_change(os.path.join('reports', 'figures', 'gene_sa.svg'))
    """

    def __init__(self, negative_ppi_file_path: str, stats_file_path: str):
        """
        Constructor
        """
        super().__init__(negative_ppi_file_path, stats_file_path)

    def _feature_init(self):
        """
        Create features for plotting surface area
        """
        self._title = "Boxplot of Change in Surface Area"
        self._ylabel = "$\Delta$ Surface Area ($\AA$)"  # noqa: W605

        # Create data object
        stats_df = pd.read_csv(
            self._stats_file_path, usecols=["symbol", "avg_sa_i_1", "avg_sa_i_2"]
        )
        neg_sa_i_1 = []
        neg_sa_i_2 = []
        pos_sa_i_1 = []
        pos_sa_i_2 = []
        for _, row in stats_df.iterrows():
            sym = row["symbol"]
            avg_sa_i_1 = row["avg_sa_i_1"]
            avg_sa_i_2 = row["avg_sa_i_2"]
            if sym in self._negative_ppis:
                neg_sa_i_1.append(avg_sa_i_1)
                neg_sa_i_2.append(avg_sa_i_2)
            else:
                pos_sa_i_1.append(avg_sa_i_1)
                pos_sa_i_2.append(avg_sa_i_2)
        data = {
            self._ylabel: pos_sa_i_1 + neg_sa_i_1 + pos_sa_i_2 + neg_sa_i_2,
            "Class": ["Cancer Driver (+)"] * len(pos_sa_i_1)
            + ["Cancer Driver (-)"] * len(neg_sa_i_1)
            + ["Partner (+)"] * len(pos_sa_i_2)
            + ["Partner (-)"] * len(neg_sa_i_2),
        }
        self._data = pd.DataFrame(data)
        self._pairs = [
            ("Cancer Driver (+)", "Cancer Driver (-)"),
            ("Partner (+)", "Partner (-)"),
        ]
        self._colour_palette = {
            "Cancer Driver (+)": "#36F1CD",
            "Cancer Driver (-)": "#D770C7",
            "Partner (+)": "#36F1CD",
            "Partner (-)": "#D770C7",
        }

    def plot_biggest_sa_change(
        self,
        save_path: str,
        top_n: int = 10,
        interaction_threshold: int = 10,
        figsize: tuple[int, int] = (8, 6),
        color="#42a5f5",
    ):
        """
        Plot the top n genes with the greatest change in surface area from
        positive to negative cases.

        Parameters
        ----------
        save_path : str
            The path to save the output graph
        top_n : int
            The number of genes you want to display (e.g., top_n = 6 means you want to display the top 6 genes).
            Defaults to 10.
        interaction_threshold : int
            The number of cases you want to use as the threshold for a valid gene observation. Lower threshold
            numbers might lead to more drastic differences between positive and negative cases but these large
            differences will be subject to small sample size bias (almost certainly). Defaults to 10.
        figsize : tuple[int, int]
            The size of the figure you wish to plot. Defaults to (8, 6).
        color : str
            The color of the bars in the bar graph
        """
        stats_df = pd.read_csv(self._stats_file_path, usecols=["symbol", "avg_sa_i_1"])
        pos_sa: defaultdict[str, list[float]] = defaultdict(list[float])
        neg_sa: defaultdict[str, list[float]] = defaultdict(list[float])
        interaction_counts: defaultdict[str, int] = defaultdict(int)

        for _, row in stats_df.iterrows():
            sym = row["symbol"]
            driver = sym.split("_")[0]
            sa = row["avg_sa_i_1"]
            interaction_counts[driver] += 1
            if sym in self._negative_ppis:
                neg_sa[sym.split("_")[0]].append(sa)
            else:
                pos_sa[sym.split("_")[0]].append(sa)

        avg_pos_sa = {
            gene: sum(sas) / len(sas)
            for gene, sas in pos_sa.items()
            if interaction_counts[gene] > interaction_threshold
        }
        avg_neg_sa = {
            gene: sum(sas) / len(sas)
            for gene, sas in neg_sa.items()
            if interaction_counts[gene] > interaction_threshold
        }
        sa_changes = {
            gene: abs(avg_pos_sa[gene] - avg_neg_sa.get(gene, 0))
            for gene in avg_pos_sa.keys()
        }
        top_genes = sorted(
            sa_changes, key=lambda gene: abs(sa_changes[gene]), reverse=True
        )[:top_n]
        genes = top_genes
        changes = [sa_changes[gene] for gene in top_genes]

        # Plot figure
        plt.figure(figsize=figsize)
        plt.bar(genes, changes, color=color)
        plt.xlabel("Gene")
        plt.ylabel("|$\Delta$ SA| Between Positive and Negative Cases")  # noqa: W605
        plt.title(f"Top {top_n} Cancer Drivers with Greatest Change in Surface Area")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(save_path, transparent=True)
