"""
Heatmap visualization module for mass imaging analysis.

This module provides functions for creating publication-ready heatmaps with:
- Data normalization/standardization
- Hierarchical clustering (rows and columns)
- Dendrogram visualization
- Customizable color schemes
- Data export for external tools
"""

import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from pathlib import Path
from typing import Optional, Tuple, Dict, List


class HeatmapPlotter:
    """
    A class for creating customizable heatmaps with clustering and normalization.
    """
    
    def __init__(self, data: pd.DataFrame):
        """
        Initialize the HeatmapPlotter.
        
        Args:
            data: DataFrame with m/z bins as rows and samples as columns
        """
        self.data = data.copy()
        self.normalized_data = None
        self.row_linkage = None
        self.col_linkage = None
        
    def normalize_data(
        self,
        method: str = 'zscore',
        axis: int = 0,
        log_transform: bool = False
    ) -> pd.DataFrame:
        """
        Normalize or standardize the data.
        
        Args:
            method: Normalization method
                - 'zscore': Z-score normalization (mean=0, std=1)
                - 'minmax': Min-max scaling (0 to 1)
                - 'robust': Robust scaling (median and IQR)
                - 'none': No normalization
            axis: Axis along which to normalize (0=rows, 1=columns)
            log_transform: Apply log2 transformation before normalization
        
        Returns:
            Normalized DataFrame
        """
        data = self.data.copy()
        
        # Apply log transformation if requested
        if log_transform:
            # Add small constant to avoid log(0)
            data = np.log2(data + 1)
        
        if method == 'none':
            self.normalized_data = data
            return data
        
        # Select scaler
        if method == 'zscore':
            scaler = StandardScaler()
        elif method == 'minmax':
            scaler = MinMaxScaler()
        elif method == 'robust':
            scaler = RobustScaler()
        else:
            raise ValueError(f"Unknown normalization method: {method}")
        
        # Apply normalization
        if axis == 0:  # Normalize each row (across samples)
            normalized = scaler.fit_transform(data.T).T
        else:  # Normalize each column (across features)
            normalized = scaler.fit_transform(data)
        
        self.normalized_data = pd.DataFrame(
            normalized,
            index=data.index,
            columns=data.columns
        )
        
        return self.normalized_data
    
    def cluster_data(
        self,
        cluster_rows: bool = True,
        cluster_cols: bool = True,
        method: str = 'average',
        metric: str = 'euclidean'
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Perform hierarchical clustering on rows and/or columns.
        
        Args:
            cluster_rows: Whether to cluster rows
            cluster_cols: Whether to cluster columns
            method: Linkage method ('average', 'complete', 'single', 'ward')
            metric: Distance metric ('euclidean', 'correlation', 'cosine', etc.)
        
        Returns:
            Tuple of (row_linkage, col_linkage)
        """
        data = self.normalized_data if self.normalized_data is not None else self.data
        
        # Cluster rows
        if cluster_rows:
            self.row_linkage = linkage(
                pdist(data, metric=metric),
                method=method
            )
        else:
            self.row_linkage = None
        
        # Cluster columns
        if cluster_cols:
            self.col_linkage = linkage(
                pdist(data.T, metric=metric),
                method=method
            )
        else:
            self.col_linkage = None
        
        return self.row_linkage, self.col_linkage
    
    def plot_heatmap(
        self,
        figsize: Tuple[int, int] = (12, 10),
        cmap: str = 'RdBu_r',
        center: Optional[float] = 0,
        show_row_dendrogram: bool = True,
        show_col_dendrogram: bool = True,
        dendrogram_ratio: float = 0.15,
        cbar_pos: Tuple[float, float, float, float] = (0.92, 0.2, 0.02, 0.6),
        row_cluster: bool = True,
        col_cluster: bool = True,
        cluster_method: str = 'average',
        cluster_metric: str = 'euclidean',
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        xlabel: str = 'Samples',
        ylabel: str = 'm/z bins',
        title: Optional[str] = None,
        annot: bool = False,
        fmt: str = '.2f',
        linewidths: float = 0,
        linecolor: str = 'white',
        xticklabels: bool = True,
        yticklabels: bool = True,
        cbar_label: str = 'Normalized Intensity'
    ) -> plt.Figure:
        """
        Create a clustered heatmap with dendrograms.
        
        Args:
            figsize: Figure size (width, height)
            cmap: Colormap name
            center: Value to center the colormap at
            show_row_dendrogram: Show row dendrogram
            show_col_dendrogram: Show column dendrogram
            dendrogram_ratio: Size ratio of dendrogram to heatmap
            cbar_pos: Colorbar position (left, bottom, width, height)
            row_cluster: Cluster rows
            col_cluster: Cluster columns
            cluster_method: Clustering linkage method
            cluster_metric: Distance metric for clustering
            vmin: Minimum value for colormap
            vmax: Maximum value for colormap
            xlabel: X-axis label
            ylabel: Y-axis label
            title: Plot title
            annot: Annotate cells with values
            fmt: Format string for annotations
            linewidths: Width of lines between cells
            linecolor: Color of lines between cells
            xticklabels: Show x-axis tick labels
            yticklabels: Show y-axis tick labels
            cbar_label: Colorbar label
        
        Returns:
            matplotlib Figure object
        """
        # Perform clustering if requested
        if row_cluster or col_cluster:
            self.cluster_data(
                cluster_rows=row_cluster,
                cluster_cols=col_cluster,
                method=cluster_method,
                metric=cluster_metric
            )
        
        # Get data to plot
        data = self.normalized_data if self.normalized_data is not None else self.data
        
        # Reorder data based on clustering
        if row_cluster and self.row_linkage is not None:
            row_order = leaves_list(self.row_linkage)
            data = data.iloc[row_order, :]
        
        if col_cluster and self.col_linkage is not None:
            col_order = leaves_list(self.col_linkage)
            data = data.iloc[:, col_order]
        
        # Create figure with subplots for dendrograms and heatmap
        fig = plt.figure(figsize=figsize)
        
        # Calculate subplot positions
        dendrogram_height = dendrogram_ratio if show_row_dendrogram else 0
        dendrogram_width = dendrogram_ratio if show_col_dendrogram else 0
        
        # Create grid spec
        from matplotlib.gridspec import GridSpec
        
        if show_row_dendrogram and show_col_dendrogram:
            gs = GridSpec(
                2, 2,
                width_ratios=[dendrogram_width, 1],
                height_ratios=[dendrogram_height, 1],
                hspace=0.02, wspace=0.02
            )
            ax_col_dendr = fig.add_subplot(gs[0, 1])
            ax_row_dendr = fig.add_subplot(gs[1, 0])
            ax_heatmap = fig.add_subplot(gs[1, 1])
        elif show_row_dendrogram:
            gs = GridSpec(
                1, 2,
                width_ratios=[dendrogram_width, 1],
                wspace=0.02
            )
            ax_row_dendr = fig.add_subplot(gs[0, 0])
            ax_heatmap = fig.add_subplot(gs[0, 1])
            ax_col_dendr = None
        elif show_col_dendrogram:
            gs = GridSpec(
                2, 1,
                height_ratios=[dendrogram_height, 1],
                hspace=0.02
            )
            ax_col_dendr = fig.add_subplot(gs[0, 0])
            ax_heatmap = fig.add_subplot(gs[1, 0])
            ax_row_dendr = None
        else:
            ax_heatmap = fig.add_subplot(111)
            ax_row_dendr = None
            ax_col_dendr = None
        
        # Plot column dendrogram
        if show_col_dendrogram and col_cluster and self.col_linkage is not None:
            dendrogram(
                self.col_linkage,
                ax=ax_col_dendr,
                orientation='top',
                no_labels=True,
                color_threshold=0
            )
            ax_col_dendr.axis('off')
        
        # Plot row dendrogram
        if show_row_dendrogram and row_cluster and self.row_linkage is not None:
            dendrogram(
                self.row_linkage,
                ax=ax_row_dendr,
                orientation='left',
                no_labels=True,
                color_threshold=0
            )
            ax_row_dendr.axis('off')
        
        # Plot heatmap
        sns.heatmap(
            data,
            ax=ax_heatmap,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            annot=annot,
            fmt=fmt,
            linewidths=linewidths,
            linecolor=linecolor,
            xticklabels=xticklabels,
            yticklabels=yticklabels,
            cbar_kws={'label': cbar_label},
            cbar_ax=None  # Let seaborn handle colorbar
        )
        
        # Set labels
        ax_heatmap.set_xlabel(xlabel, fontweight='bold', fontsize=12)
        ax_heatmap.set_ylabel(ylabel, fontweight='bold', fontsize=12)
        
        # Rotate labels
        ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=45, ha='right')
        ax_heatmap.set_yticklabels(ax_heatmap.get_yticklabels(), rotation=0)
        
        # Set title
        if title:
            fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
        
        plt.tight_layout()
        
        return fig
    
    def get_clustered_data(self) -> pd.DataFrame:
        """
        Get the data in clustered order.
        
        Returns:
            DataFrame with rows and columns reordered by clustering
        """
        data = self.normalized_data if self.normalized_data is not None else self.data
        
        if self.row_linkage is not None:
            row_order = leaves_list(self.row_linkage)
            data = data.iloc[row_order, :]
        
        if self.col_linkage is not None:
            col_order = leaves_list(self.col_linkage)
            data = data.iloc[:, col_order]
        
        return data
    
    def export_data(
        self,
        output_path: str,
        include_original: bool = True,
        include_normalized: bool = True,
        include_clustered: bool = True
    ) -> None:
        """
        Export data to CSV file(s).
        
        Args:
            output_path: Base path for output files
            include_original: Export original data
            include_normalized: Export normalized data
            include_clustered: Export clustered data
        """
        output_path = Path(output_path)
        base_name = output_path.stem
        output_dir = output_path.parent
        
        if include_original:
            original_path = output_dir / f"{base_name}_original.csv"
            self.data.to_csv(original_path, float_format='%.4f')
            print(f"✓ Saved original data: {original_path}")
        
        if include_normalized and self.normalized_data is not None:
            normalized_path = output_dir / f"{base_name}_normalized.csv"
            self.normalized_data.to_csv(normalized_path, float_format='%.4f')
            print(f"✓ Saved normalized data: {normalized_path}")
        
        if include_clustered:
            clustered_data = self.get_clustered_data()
            clustered_path = output_dir / f"{base_name}_clustered.csv"
            clustered_data.to_csv(clustered_path, float_format='%.4f')
            print(f"✓ Saved clustered data: {clustered_path}")


def load_integrated_results(
    file_path: str,
    sample_columns: Optional[List[str]] = None,
    index_column: str = 'm_z_bin'
) -> pd.DataFrame:
    """
    Load integrated results CSV and extract sample intensity data.
    
    Args:
        file_path: Path to integrated_results CSV file
        sample_columns: List of sample column names to extract
                       If None, auto-detect columns matching pattern {group}_{n}
        index_column: Column to use as row index (default: 'm_z_bin')
    
    Returns:
        DataFrame with m/z bins as rows and samples as columns
    """
    df = pd.read_csv(file_path)
    
    # Auto-detect sample columns if not provided
    if sample_columns is None:
        # Find columns matching pattern: group_number (e.g., saline_1, glyoxylate_2)
        sample_columns = [col for col in df.columns 
                         if '_' in col and col.split('_')[-1].isdigit()
                         and not col.startswith('log2FC')
                         and not col.startswith('FC_')
                         and not col.endswith('_mean')
                         and not col.endswith('_sd')]
    
    # Extract sample data
    if index_column in df.columns:
        data = df.set_index(index_column)[sample_columns]
    else:
        data = df[sample_columns]
    
    return data


def create_heatmap_from_file(
    file_path: str,
    output_path: Optional[str] = None,
    normalization: str = 'zscore',
    log_transform: bool = False,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    cluster_method: str = 'average',
    cluster_metric: str = 'euclidean',
    show_row_dendrogram: bool = True,
    show_col_dendrogram: bool = True,
    figsize: Tuple[int, int] = (12, 10),
    cmap: str = 'RdBu_r',
    title: Optional[str] = None,
    **kwargs
) -> Tuple[plt.Figure, HeatmapPlotter]:
    """
    Convenience function to create a heatmap from an integrated results file.
    
    Args:
        file_path: Path to integrated_results CSV file
        output_path: Path to save the figure (optional)
        normalization: Normalization method ('zscore', 'minmax', 'robust', 'none')
        log_transform: Apply log2 transformation before normalization
        cluster_rows: Cluster rows (m/z bins)
        cluster_cols: Cluster columns (samples)
        cluster_method: Clustering linkage method
        cluster_metric: Distance metric for clustering
        show_row_dendrogram: Show row dendrogram
        show_col_dendrogram: Show column dendrogram
        figsize: Figure size
        cmap: Colormap
        title: Plot title
        **kwargs: Additional arguments passed to plot_heatmap
    
    Returns:
        Tuple of (figure, HeatmapPlotter instance)
    """
    # Load data
    data = load_integrated_results(file_path)
    
    # Create plotter
    plotter = HeatmapPlotter(data)
    
    # Normalize
    plotter.normalize_data(method=normalization, axis=0, log_transform=log_transform)
    
    # Create plot
    fig = plotter.plot_heatmap(
        figsize=figsize,
        cmap=cmap,
        title=title,
        show_row_dendrogram=show_row_dendrogram,
        show_col_dendrogram=show_col_dendrogram,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        cluster_method=cluster_method,
        cluster_metric=cluster_metric,
        **kwargs
    )
    
    # Save if output path provided
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"✓ Saved heatmap: {output_path}")
    
    return fig, plotter
