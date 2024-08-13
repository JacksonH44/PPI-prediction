"""
A class for plotting changes in features of AF-M folded complexes
from the positive and negative set in the dataset.
"""

from abc import ABC, abstractmethod

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator

class Plotter(ABC):
    """
    An abstract class for plotting the differences in features
    stratified by split (positive or negative).

    Attributes
    ----------
    _negative_ppi_file_path : str
        The path to the file that holds all negative PPIs
    _stats_file_path : str
        The path to the file that holds all stats from the complexes
    _data : pd.DataFrame
        The dataframe representing the data that will be plotted
    _ylabel : str
        The label of the y-axis of the boxplot
    _title : str
        The title of the boxplot
    _pairs : list[tuple[str]]
        Pairs of categories of data (e.g., ('Cancer Driver (+)', 'Cancer Driver (-)'))
        to be compared on a statistical significance scale with a t-test.
    _columns : list[str]
        Columns from the ColabFold stats file to extract when creating _data
    _colour_palette : dict[str, str]
        A map of category, colour pairs for boxplot plotting
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
        self._negative_ppis = self._create_negative_ppis()
        self._data : pd.DataFrame = None
        self._ylabel : str = None
        self._title : str = None
        self._pairs : list[tuple[str]] = None
        self._columns : list[str] = None
        self._colour_palette : dict[str, str] = None
        self._feature_init()

    @abstractmethod
    def _feature_init(self):
        """
        Should initialize self._data, self._ylabel, self._pairs,
        self._title, self._columns, and self._colour_palette.
        """

    def _create_negative_ppis(self) -> set[str]:
        """
        Generate a set of negative protein-protein interactions
        used to split the dataset when plotting.

        Returns
        -------
        negative_ppis : set[str]
            A set of all negative PPIs in the dataset
        """
        negative_dataset = pd.read_csv(self._negative_ppi_file_path)
        negative_ppis = set()
        for _, row in negative_dataset.iterrows():
            s = f"{row['gene_symbol_a']}_{row['gene_symbol_b']}"
            negative_ppis.add(s)
        return negative_ppis
    
    def plot_feature(self, save_path: str, figsize: tuple[int, int] = (8, 6)):
        """
        Create a boxplot with a stripplot overlaid of the difference in a
        feature from a positive to negative case.

        Parameters
        ----------
        save_path : str
            Path to which you wish to save the plot to.
        figsize : tuple[int, int]
            A tuple of (width, height) that specifies how large the figure
            should be. Defaults to (8, 6).
        """
        plt.figure(figsize=figsize)

        # Create the boxplot
        ax = sns.stripplot(x='Class', y=self._ylabel, data=self._data, palette=self._colour_palette, alpha=0.4)
        sns.boxplot(x='Class', y=self._ylabel, data=self._data, showfliers=False, ax=ax, palette=self._colour_palette)
        plt.title(self._title)
        plt.ylabel(self._ylabel)

        # Annotate the median values
        median_values = self._data.groupby('Class')[self._ylabel].median().to_dict()
        for i, classification in enumerate(self._data['Class'].unique()):
            median_val = median_values[classification]
            ax.annotate(f'{median_val:.2f}', 
                    xy=(i, median_val), 
                    xycoords='data',
                    xytext=(0, 10), 
                    textcoords='offset points',
                    ha='center', va='center', 
                    fontsize=10, color='black',
                    bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))
        
        # Annotate the statistical significance
        annotater = Annotator(ax, self._pairs, data=self._data, x="Class", y=self._ylabel)
        annotater.configure(test="t-test_ind", text_format="star")
        annotater.apply_and_annotate()

        # Plot figure
        plt.savefig(save_path, transparent=True)
