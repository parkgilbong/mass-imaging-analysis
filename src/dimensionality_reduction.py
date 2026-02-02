"""
Dimensionality reduction module for mass imaging analysis.

This module provides functions for dimensionality reduction and visualization:
- PCA (Principal Component Analysis)
- t-SNE (t-Distributed Stochastic Neighbor Embedding)
- UMAP (Uniform Manifold Approximation and Projection)
- Data normalization/standardization
- 2D scatter plot visualization
- Data export for external tools
"""

import warnings
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Union


class DimensionalityReducer:
    """
    A class for performing dimensionality reduction and visualization.
    """
    
    def __init__(self, data: pd.DataFrame, sample_groups: Optional[Dict[str, str]] = None):
        """
        Initialize the DimensionalityReducer.
        
        Args:
            data: DataFrame with samples as rows and features (m/z bins) as columns
            sample_groups: Dictionary mapping sample names to group labels
        """
        self.data = data.copy()
        self.normalized_data = None
        self.reduced_data = None
        self.algorithm = None
        self.model = None
        
        # Auto-detect groups from sample names if not provided
        if sample_groups is None:
            self.sample_groups = self._auto_detect_groups()
        else:
            self.sample_groups = sample_groups
    
    def _auto_detect_groups(self) -> Dict[str, str]:
        """
        Auto-detect group labels from sample names.
        Assumes sample names follow pattern: {group}_{number}
        
        Returns:
            Dictionary mapping sample names to group labels
        """
        groups = {}
        for sample in self.data.index:
            # Extract group name (everything before last underscore)
            if '_' in sample:
                group = sample.rsplit('_', 1)[0]
            else:
                group = 'unknown'
            groups[sample] = group
        return groups
    
    def normalize_data(
        self,
        method: str = 'zscore',
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
            log_transform: Apply log2 transformation before normalization
        
        Returns:
            Normalized DataFrame
        """
        data = self.data.copy()
        
        # Apply log transformation if requested
        if log_transform:
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
        
        # Apply normalization (normalize across features for each sample)
        normalized = scaler.fit_transform(data)
        
        self.normalized_data = pd.DataFrame(
            normalized,
            index=data.index,
            columns=data.columns
        )
        
        return self.normalized_data
    
    def apply_pca(
        self,
        n_components: int = 2,
        random_state: int = 42
    ) -> pd.DataFrame:
        """
        Apply PCA (Principal Component Analysis).
        
        Args:
            n_components: Number of principal components
            random_state: Random seed for reproducibility
        
        Returns:
            DataFrame with reduced dimensions
        """
        data = self.normalized_data if self.normalized_data is not None else self.data
        
        # Apply PCA
        pca = PCA(n_components=n_components, random_state=random_state)
        reduced = pca.fit_transform(data)
        
        # Store model and results
        self.model = pca
        self.algorithm = 'PCA'
        
        # Create DataFrame
        columns = [f'PC{i+1}' for i in range(n_components)]
        self.reduced_data = pd.DataFrame(
            reduced,
            index=data.index,
            columns=columns
        )
        
        return self.reduced_data
    
    def apply_tsne(
        self,
        n_components: int = 2,
        perplexity: float = 30.0,
        learning_rate: Union[float, str] = 'auto',
        n_iter: int = 1000,
        random_state: int = 42
    ) -> pd.DataFrame:
        """
        Apply t-SNE (t-Distributed Stochastic Neighbor Embedding).
        
        Args:
            n_components: Number of dimensions (usually 2)
            perplexity: Perplexity parameter (5-50 recommended)
            learning_rate: Learning rate ('auto' or float, 10-1000)
            n_iter: Number of iterations (250-1000 recommended)
            random_state: Random seed for reproducibility
        
        Returns:
            DataFrame with reduced dimensions
        """
        data = self.normalized_data if self.normalized_data is not None else self.data
        
        # Apply t-SNE
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            n_iter=n_iter,
            random_state=random_state
        )
        reduced = tsne.fit_transform(data)
        
        # Store model and results
        self.model = tsne
        self.algorithm = 't-SNE'
        
        # Create DataFrame
        columns = [f't-SNE{i+1}' for i in range(n_components)]
        self.reduced_data = pd.DataFrame(
            reduced,
            index=data.index,
            columns=columns
        )
        
        return self.reduced_data
    
    def apply_umap(
        self,
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        metric: str = 'euclidean',
        random_state: int = 42
    ) -> pd.DataFrame:
        """
        Apply UMAP (Uniform Manifold Approximation and Projection).
        
        Args:
            n_components: Number of dimensions (usually 2)
            n_neighbors: Number of neighbors (2-100, default 15)
            min_dist: Minimum distance (0.0-0.99, default 0.1)
            metric: Distance metric ('euclidean', 'manhattan', 'cosine', etc.)
            random_state: Random seed for reproducibility
        
        Returns:
            DataFrame with reduced dimensions
        """
        try:
            import umap
        except ImportError:
            raise ImportError(
                "UMAP is not installed. Install it with: pip install umap-learn"
            )
        
        data = self.normalized_data if self.normalized_data is not None else self.data
        
        # Apply UMAP
        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            metric=metric,
            random_state=random_state
        )
        reduced = reducer.fit_transform(data)
        
        # Store model and results
        self.model = reducer
        self.algorithm = 'UMAP'
        
        # Create DataFrame
        columns = [f'UMAP{i+1}' for i in range(n_components)]
        self.reduced_data = pd.DataFrame(
            reduced,
            index=data.index,
            columns=columns
        )
        
        return self.reduced_data
    
    def plot_2d(
        self,
        figsize: Tuple[int, int] = (10, 8),
        point_size: int = 100,
        alpha: float = 0.7,
        show_labels: bool = True,
        label_fontsize: int = 9,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        legend_loc: str = 'best',
        palette: Optional[str] = 'Set2',
        show_variance: bool = True
    ) -> plt.Figure:
        """
        Create a 2D scatter plot of the reduced data.
        
        Args:
            figsize: Figure size (width, height)
            point_size: Size of scatter points
            alpha: Transparency of points
            show_labels: Show sample labels
            label_fontsize: Font size for labels
            title: Plot title (auto-generated if None)
            xlabel: X-axis label (auto-generated if None)
            ylabel: Y-axis label (auto-generated if None)
            legend_loc: Legend location
            palette: Color palette name
            show_variance: Show variance explained (PCA only)
        
        Returns:
            matplotlib Figure object
        """
        if self.reduced_data is None:
            raise ValueError("No reduced data available. Run a dimensionality reduction method first.")
        
        if self.reduced_data.shape[1] < 2:
            raise ValueError("Need at least 2 dimensions for 2D plot")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get unique groups and assign colors
        unique_groups = sorted(set(self.sample_groups.values()))
        colors = sns.color_palette(palette, n_colors=len(unique_groups))
        group_colors = dict(zip(unique_groups, colors))
        
        # Plot each group
        for group in unique_groups:
            # Get samples in this group
            group_samples = [s for s, g in self.sample_groups.items() if g == group]
            group_data = self.reduced_data.loc[group_samples]
            
            # Plot
            ax.scatter(
                group_data.iloc[:, 0],
                group_data.iloc[:, 1],
                c=[group_colors[group]],
                s=point_size,
                alpha=alpha,
                label=group,
                edgecolors='black',
                linewidths=0.5
            )
            
            # Add labels if requested
            if show_labels:
                for idx, row in group_data.iterrows():
                    ax.annotate(
                        idx,
                        (row.iloc[0], row.iloc[1]),
                        fontsize=label_fontsize,
                        alpha=0.7,
                        xytext=(5, 5),
                        textcoords='offset points'
                    )
        
        # Set labels
        if xlabel is None:
            xlabel = self.reduced_data.columns[0]
            if self.algorithm == 'PCA' and show_variance and hasattr(self.model, 'explained_variance_ratio_'):
                var1 = self.model.explained_variance_ratio_[0] * 100
                xlabel = f'{xlabel} ({var1:.1f}%)'
        
        if ylabel is None:
            ylabel = self.reduced_data.columns[1]
            if self.algorithm == 'PCA' and show_variance and hasattr(self.model, 'explained_variance_ratio_'):
                var2 = self.model.explained_variance_ratio_[1] * 100
                ylabel = f'{ylabel} ({var2:.1f}%)'
        
        ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        
        # Set title
        if title is None:
            title = f'{self.algorithm} Analysis'
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # Add legend
        ax.legend(loc=legend_loc, frameon=True, fancybox=True, shadow=True)
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle=':')
        
        plt.tight_layout()
        
        return fig
    
    def plot_variance_explained(
        self,
        n_components: int = 10,
        figsize: Tuple[int, int] = (10, 6)
    ) -> plt.Figure:
        """
        Plot variance explained by principal components (PCA only).
        
        Args:
            n_components: Number of components to show
            figsize: Figure size
        
        Returns:
            matplotlib Figure object
        """
        if self.algorithm != 'PCA':
            raise ValueError("Variance explained plot is only available for PCA")
        
        if not hasattr(self.model, 'explained_variance_ratio_'):
            raise ValueError("PCA model does not have explained_variance_ratio_")
        
        # Get variance explained
        var_ratio = self.model.explained_variance_ratio_
        n_show = min(n_components, len(var_ratio))
        
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        
        # Individual variance explained
        ax1.bar(
            range(1, n_show + 1),
            var_ratio[:n_show] * 100,
            alpha=0.7,
            color='steelblue',
            edgecolor='black'
        )
        ax1.set_xlabel('Principal Component', fontweight='bold')
        ax1.set_ylabel('Variance Explained (%)', fontweight='bold')
        ax1.set_title('Individual Variance Explained', fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Cumulative variance explained
        cumsum = np.cumsum(var_ratio[:n_show]) * 100
        ax2.plot(
            range(1, n_show + 1),
            cumsum,
            marker='o',
            linewidth=2,
            markersize=8,
            color='steelblue'
        )
        ax2.axhline(y=80, color='red', linestyle='--', alpha=0.5, label='80%')
        ax2.axhline(y=90, color='orange', linestyle='--', alpha=0.5, label='90%')
        ax2.set_xlabel('Principal Component', fontweight='bold')
        ax2.set_ylabel('Cumulative Variance Explained (%)', fontweight='bold')
        ax2.set_title('Cumulative Variance Explained', fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        
        return fig
    
    def export_data(
        self,
        output_path: str,
        include_original: bool = True,
        include_normalized: bool = True,
        include_reduced: bool = True
    ) -> None:
        """
        Export data to CSV file(s).
        
        Args:
            output_path: Base path for output files
            include_original: Export original data
            include_normalized: Export normalized data
            include_reduced: Export reduced data
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
        
        if include_reduced and self.reduced_data is not None:
            reduced_path = output_dir / f"{base_name}_{self.algorithm}_reduced.csv"
            
            # Add group information
            reduced_with_groups = self.reduced_data.copy()
            reduced_with_groups.insert(0, 'group', [self.sample_groups[s] for s in reduced_with_groups.index])
            
            reduced_with_groups.to_csv(reduced_path, float_format='%.4f')
            print(f"✓ Saved {self.algorithm} reduced data: {reduced_path}")


def load_integrated_results_for_dr(
    file_path: str,
    sample_columns: Optional[List[str]] = None
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Load integrated results CSV and prepare for dimensionality reduction.
    
    Args:
        file_path: Path to integrated_results CSV file
        sample_columns: List of sample column names to extract
                       If None, auto-detect columns matching pattern {group}_{n}
    
    Returns:
        Tuple of (data DataFrame, sample_groups dictionary)
        - data: Samples as rows, m/z bins as columns
        - sample_groups: Dictionary mapping sample names to group labels
    """
    df = pd.read_csv(file_path)
    
    # Auto-detect sample columns if not provided
    if sample_columns is None:
        sample_columns = [col for col in df.columns 
                         if '_' in col and col.split('_')[-1].isdigit()
                         and not col.startswith('log2FC')
                         and not col.startswith('FC_')
                         and not col.endswith('_mean')
                         and not col.endswith('_sd')]
    
    # Extract sample data (transpose so samples are rows)
    data = df[sample_columns].T
    data.columns = df['m_z_bin'].astype(str)
    
    # Create sample groups dictionary
    sample_groups = {}
    for sample in data.index:
        if '_' in sample:
            group = sample.rsplit('_', 1)[0]
        else:
            group = 'unknown'
        sample_groups[sample] = group
    
    return data, sample_groups


def quick_pca(
    file_path: str,
    output_path: Optional[str] = None,
    normalization: str = 'zscore',
    log_transform: bool = False,
    n_components: int = 2,
    **kwargs
) -> Tuple[plt.Figure, DimensionalityReducer]:
    """
    Convenience function to quickly create a PCA plot.
    
    Args:
        file_path: Path to integrated_results CSV file
        output_path: Path to save the figure (optional)
        normalization: Normalization method
        log_transform: Apply log transformation
        n_components: Number of principal components
        **kwargs: Additional arguments passed to plot_2d
    
    Returns:
        Tuple of (figure, DimensionalityReducer instance)
    """
    # Load data
    data, sample_groups = load_integrated_results_for_dr(file_path)
    
    # Create reducer
    reducer = DimensionalityReducer(data, sample_groups)
    
    # Normalize
    reducer.normalize_data(method=normalization, log_transform=log_transform)
    
    # Apply PCA
    reducer.apply_pca(n_components=n_components)
    
    # Plot
    fig = reducer.plot_2d(**kwargs)
    
    # Save if output path provided
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"✓ Saved PCA plot: {output_path}")
    
    return fig, reducer
