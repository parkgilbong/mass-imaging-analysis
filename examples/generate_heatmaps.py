"""
Example script demonstrating heatmap generation.

This script shows how to use the heatmap_plot module programmatically.
"""

import sys
from pathlib import Path

# Add src directory to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from heatmap_plot import create_heatmap_from_file

def main():
    """Generate example heatmaps with different settings."""
    
    # Input file
    input_file = 'output/4_groups/integrated_results_cortex.csv'
    
    # Output directory
    output_dir = Path('output/4_groups/heatmaps')
    output_dir.mkdir(exist_ok=True, parents=True)
    
    print("="*60)
    print("Heatmap Generation Examples")
    print("="*60)
    
    # Example 1: Basic heatmap with default settings
    print("\n1. Creating basic heatmap...")
    fig1, plotter1 = create_heatmap_from_file(
        file_path=input_file,
        output_path=str(output_dir / 'heatmap_basic.png'),
        normalization='zscore',
        cluster_rows=True,
        cluster_cols=True,
        title='Basic Heatmap (Z-score normalized)'
    )
    
    # Export data
    plotter1.export_data(
        str(output_dir / 'heatmap_basic_data'),
        include_original=True,
        include_normalized=True,
        include_clustered=True
    )
    
    # Example 2: Heatmap with log transformation
    print("\n2. Creating heatmap with log transformation...")
    fig2, plotter2 = create_heatmap_from_file(
        file_path=input_file,
        output_path=str(output_dir / 'heatmap_log.png'),
        normalization='zscore',
        log_transform=True,
        cluster_rows=True,
        cluster_cols=True,
        title='Heatmap (Log2 + Z-score)'
    )
    
    # Example 3: Heatmap with different colormap
    print("\n3. Creating heatmap with viridis colormap...")
    fig3, plotter3 = create_heatmap_from_file(
        file_path=input_file,
        output_path=str(output_dir / 'heatmap_viridis.png'),
        normalization='minmax',
        cluster_rows=True,
        cluster_cols=True,
        cmap='viridis',
        center=None,  # No centering for sequential colormap
        title='Heatmap (Min-Max + Viridis)',
        cbar_label='Normalized Intensity (0-1)'
    )
    
    # Example 4: Heatmap without column clustering (preserve sample order)
    print("\n4. Creating heatmap without column clustering...")
    fig4, plotter4 = create_heatmap_from_file(
        file_path=input_file,
        output_path=str(output_dir / 'heatmap_no_col_cluster.png'),
        normalization='zscore',
        cluster_rows=True,
        cluster_cols=False,  # Don't cluster columns
        show_col_dendrogram=False,
        title='Heatmap (Rows clustered, columns ordered)'
    )
    
    # Example 5: Heatmap with custom clustering
    print("\n5. Creating heatmap with ward clustering...")
    fig5, plotter5 = create_heatmap_from_file(
        file_path=input_file,
        output_path=str(output_dir / 'heatmap_ward.png'),
        normalization='zscore',
        cluster_rows=True,
        cluster_cols=True,
        cluster_method='ward',  # Ward linkage
        cluster_metric='euclidean',
        title='Heatmap (Ward clustering)'
    )
    
    print("\n" + "="*60)
    print("All heatmaps generated successfully!")
    print(f"Output directory: {output_dir}")
    print("="*60)

if __name__ == '__main__':
    main()
