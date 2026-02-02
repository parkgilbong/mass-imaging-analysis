"""
Example script demonstrating dimensionality reduction analysis.

This script shows how to use the dimensionality_reduction module programmatically.
"""

import sys
from pathlib import Path

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from dimensionality_reduction import (
    DimensionalityReducer,
    load_integrated_results_for_dr,
    quick_pca
)
import matplotlib.pyplot as plt

def main():
    """Generate example dimensionality reduction plots."""
    
    # Input file
    input_file = 'output/4_groups/integrated_results_cortex.csv'
    
    # Output directory
    output_dir = Path('output/4_groups/dimensionality_reduction')
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print("="*60)
    print("Dimensionality Reduction Examples")
    print("="*60)
    
    # Load data
    print("\nLoading data...")
    data, sample_groups = load_integrated_results_for_dr(input_file)
    print(f"✓ Loaded {data.shape[0]} samples × {data.shape[1]} features")
    
    # Example 1: Quick PCA
    print("\n1. Quick PCA analysis...")
    fig_pca, reducer_pca = quick_pca(
        file_path=input_file,
        output_path=str(output_dir / 'pca_quick.png'),
        normalization='zscore',
        log_transform=False,
        n_components=2,
        title='PCA Analysis (Quick)'
    )
    plt.close(fig_pca)
    
    # Export PCA data
    reducer_pca.export_data(
        str(output_dir / 'pca_data'),
        include_original=True,
        include_normalized=True,
        include_reduced=True
    )
    
    # Example 2: PCA with variance explained
    print("\n2. PCA with variance analysis...")
    reducer = DimensionalityReducer(data, sample_groups)
    reducer.normalize_data(method='zscore', log_transform=False)
    reducer.apply_pca(n_components=10)
    
    # Plot variance explained
    fig_var = reducer.plot_variance_explained(n_components=10, figsize=(12, 5))
    fig_var.savefig(output_dir / 'pca_variance.png', dpi=300, bbox_inches='tight')
    plt.close(fig_var)
    print(f"✓ Saved variance plot")
    
    # Example 3: t-SNE with different perplexities
    print("\n3. t-SNE with different perplexities...")
    perplexities = [5, 30, 50]
    
    fig_tsne, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for i, perplexity in enumerate(perplexities):
        print(f"   - Perplexity {perplexity}...")
        
        reducer_tsne = DimensionalityReducer(data, sample_groups)
        reducer_tsne.normalize_data(method='zscore', log_transform=False)
        reducer_tsne.apply_tsne(
            n_components=2,
            perplexity=perplexity,
            n_iter=1000,
            random_state=42
        )
        
        # Plot on subplot
        ax = axes[i]
        
        # Get unique groups and colors
        import seaborn as sns
        unique_groups = sorted(set(sample_groups.values()))
        colors = sns.color_palette('Set2', n_colors=len(unique_groups))
        group_colors = dict(zip(unique_groups, colors))
        
        for group in unique_groups:
            group_samples = [s for s, g in sample_groups.items() if g == group]
            group_data = reducer_tsne.reduced_data.loc[group_samples]
            
            ax.scatter(
                group_data.iloc[:, 0],
                group_data.iloc[:, 1],
                c=[group_colors[group]],
                s=100,
                alpha=0.7,
                label=group,
                edgecolors='black',
                linewidths=0.5
            )
        
        ax.set_xlabel('t-SNE1', fontweight='bold')
        ax.set_ylabel('t-SNE2', fontweight='bold')
        ax.set_title(f't-SNE (perplexity={perplexity})', fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle=':')
        
        if i == 2:
            ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    fig_tsne.savefig(output_dir / 'tsne_perplexity_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig_tsne)
    print(f"✓ Saved t-SNE comparison")
    
    # Example 4: UMAP (if available)
    print("\n4. UMAP analysis...")
    try:
        reducer_umap = DimensionalityReducer(data, sample_groups)
        reducer_umap.normalize_data(method='zscore', log_transform=False)
        reducer_umap.apply_umap(
            n_components=2,
            n_neighbors=15,
            min_dist=0.1,
            random_state=42
        )
        
        fig_umap = reducer_umap.plot_2d(
            figsize=(10, 8),
            title='UMAP Analysis',
            point_size=100,
            alpha=0.7
        )
        fig_umap.savefig(output_dir / 'umap.png', dpi=300, bbox_inches='tight')
        plt.close(fig_umap)
        print(f"✓ Saved UMAP plot")
        
    except ImportError:
        print("   ⚠ UMAP not installed (install with: pip install umap-learn)")
    
    # Example 5: Algorithm comparison
    print("\n5. Algorithm comparison...")
    fig_compare, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    algorithms = [
        ('PCA', lambda r: r.apply_pca(n_components=2, random_state=42)),
        ('t-SNE', lambda r: r.apply_tsne(perplexity=30, n_iter=1000, random_state=42)),
        ('UMAP', lambda r: r.apply_umap(n_neighbors=15, min_dist=0.1, random_state=42))
    ]
    
    for i, (algo_name, algo_func) in enumerate(algorithms):
        try:
            print(f"   - {algo_name}...")
            
            r = DimensionalityReducer(data, sample_groups)
            r.normalize_data(method='zscore', log_transform=False)
            algo_func(r)
            
            ax = axes[i]
            
            import seaborn as sns
            unique_groups = sorted(set(sample_groups.values()))
            colors = sns.color_palette('Set2', n_colors=len(unique_groups))
            group_colors = dict(zip(unique_groups, colors))
            
            for group in unique_groups:
                group_samples = [s for s, g in sample_groups.items() if g == group]
                group_data = r.reduced_data.loc[group_samples]
                
                ax.scatter(
                    group_data.iloc[:, 0],
                    group_data.iloc[:, 1],
                    c=[group_colors[group]],
                    s=100,
                    alpha=0.7,
                    label=group,
                    edgecolors='black',
                    linewidths=0.5
                )
            
            ax.set_xlabel(r.reduced_data.columns[0], fontweight='bold')
            ax.set_ylabel(r.reduced_data.columns[1], fontweight='bold')
            ax.set_title(algo_name, fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle=':')
            
            if i == 2:
                ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
                
        except ImportError:
            ax.text(0.5, 0.5, f'{algo_name}\nNot Available', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
    
    plt.tight_layout()
    fig_compare.savefig(output_dir / 'algorithm_comparison.png', dpi=300, bbox_inches='tight')
    plt.close(fig_compare)
    print(f"✓ Saved algorithm comparison")
    
    print("\n" + "="*60)
    print("All plots generated successfully!")
    print(f"Output directory: {output_dir}")
    print("="*60)

if __name__ == '__main__':
    main()
