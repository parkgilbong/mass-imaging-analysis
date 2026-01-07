import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')
warnings.filterwarnings('ignore', category=FutureWarning, module='seaborn')

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime

try:
    import parsing
    from utils.logging_utils import get_logger
except ImportError:
    try:
        from . import parsing
        from .utils.logging_utils import get_logger
    except ImportError:
        pass 

# Logger will be initialized in main() with optional log file

def plot_single_bin(ax, data_bin, m_z_bin, stats_df_bin):
    """
    Plot bar chart with individual data points and pairwise comparison brackets.
    """
    # Create bar plot with strip plot overlay
    sns.barplot(data=data_bin, x='group', y='intensity', hue='group', palette='Paired', errorbar=None, ax=ax, legend=False)
    sns.stripplot(data=data_bin, x='group', y='intensity', color='black', legend=False, s=5, ax=ax, jitter=False)
    ax.set_title(f"m/z bin: {m_z_bin}", fontsize=10)
    ax.set_xlabel('')
    ax.set_ylabel('Intensity')
    ax.tick_params(axis='x', rotation=45)
    
    # Add pairwise comparison brackets for significant pairs
    significant_pairs = stats_df_bin[stats_df_bin['significant'] == True]
    
    if not significant_pairs.empty:
        # Get group names and their x-axis positions
        groups = data_bin['group'].unique()
        group_to_x = {group: i for i, group in enumerate(groups)}
        
        # Calculate y-positions for brackets
        max_intensity = data_bin['intensity'].max()
        y_range = data_bin['intensity'].max() - data_bin['intensity'].min()
        bracket_height = y_range * 0.05  # Height of each bracket level
        base_y = max_intensity + y_range * 0.05  # Starting y position
        
        # Draw brackets for each significant pair
        for idx, (_, row) in enumerate(significant_pairs.iterrows()):
            group1 = row['group1']
            group2 = row['group2']
            
            # Get x positions for the two groups
            x1 = group_to_x.get(group1)
            x2 = group_to_x.get(group2)
            
            if x1 is not None and x2 is not None:
                # Calculate y position for this bracket (stack them vertically)
                y = base_y + idx * bracket_height * 1.5
                
                # Draw horizontal lines and vertical connectors
                ax.plot([x1, x1, x2, x2], [y, y + bracket_height*0.3, y + bracket_height*0.3, y], 
                       'k-', linewidth=1.5)
                
                # Add asterisk in the middle
                x_mid = (x1 + x2) / 2
                ax.text(x_mid, y + bracket_height*0.5, '*', 
                       ha='center', va='bottom', fontsize=14, color='red', fontweight='bold')

def generate_single_plot(data_bin, roi, m_z_bin, stats_df_bin, group_color_map, groups, output_dir):
    try:
        fig_single, ax_single = plt.subplots(figsize=(6, 5))
        plot_single_bin(ax_single, data_bin, m_z_bin, stats_df_bin)
        handles = [plt.Rectangle((0,0),1,1, color=group_color_map[group]) for group in groups]
        ax_single.legend(handles, groups, title="Groups", bbox_to_anchor=(1.05, 1), loc='upper left')
        bin_filename_safe = str(m_z_bin).replace('.', '_')
        single_plot_path = os.path.join(output_dir, f"plot_roi_{roi}_bin_{bin_filename_safe}.png")
        fig_single.savefig(single_plot_path, dpi=150, bbox_inches='tight')
        plt.close(fig_single) 
    except Exception as e:
        logger.error(f"Failed to save individual plot ({m_z_bin}): {e}")
        plt.close(fig_single)

def generate_montage_plot(df_long_roi, roi, value_vars, df_posthoc_roi, group_color_map, groups, output_dir):
    logger.info(f"Generating montage plot: {roi}")
    num_bins = len(value_vars)
    cols = int(np.ceil(np.sqrt(num_bins)))
    rows = int(np.ceil(num_bins / cols))
    
    fig = plt.figure(figsize=(cols * 5, rows * 4))
    gs = gridspec.GridSpec(rows, cols + 1, figure=fig, width_ratios=[1] * cols + [0.5])
    
    plot_axes = []
    for r in range(rows):
        for c in range(cols):
            if r * cols + c < num_bins:
                plot_axes.append(fig.add_subplot(gs[r, c]))

    for idx, (m_z_bin, ax) in enumerate(zip(value_vars, plot_axes)):
        data_bin = df_long_roi[df_long_roi['m_z_bin'] == m_z_bin]
        stats_df_bin = pd.DataFrame()
        if not df_posthoc_roi.empty:
            stats_df_bin = df_posthoc_roi[(df_posthoc_roi['m_z_bin'] == m_z_bin) & (df_posthoc_roi['roi'] == roi)]
        plot_single_bin(ax, data_bin, m_z_bin, stats_df_bin)

    legend_ax = fig.add_subplot(gs[:, -1])
    legend_ax.axis('off')
    handles = [plt.Rectangle((0,0),1,1, color=group_color_map[group]) for group in groups]
    legend_ax.legend(handles, groups, title="Groups", loc='center')

    fig.suptitle(f'MSI Intensity Analysis (ROI: {roi})', fontsize=16, y=1.02)
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    
    montage_plot_path = os.path.join(output_dir, f"plot_montage_roi_{roi}.png")
    fig.savefig(montage_plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

def generate_html_report(config, output_dir, df_main_stats, df_posthoc_stats, rois):
    logger.info("Starting HTML report generation")
    
    report_path = os.path.join(output_dir, 'analysis_report.html')
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Mass Imaging Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1, h2 {{ color: #333; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 20px; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
            .plot-container {{ margin-bottom: 40px; }}
            img {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
        </style>
    </head>
    <body>
        <h1>Mass Imaging Analysis Report</h1>
        <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <h2>Configuration</h2>
        <ul>
            <li>Output Directory: {output_dir}</li>
            <li>ROIs: {', '.join(rois)}</li>
        </ul>

        <h2>Statistical Results (Main)</h2>
        <div style="overflow-x:auto;">
            {df_main_stats.to_html(index=False, classes='table table-striped')}
        </div>
    """
    
    if not df_posthoc_stats.empty:
        html_content += f"""
        <h2>Post-hoc Results</h2>
        <div style="overflow-x:auto;">
            {df_posthoc_stats.to_html(index=False, classes='table table-striped')}
        </div>
        """
        
    html_content += "<h2>Plots</h2>"
    
    for roi in rois:
        montage_path = f"plot_montage_roi_{roi}.png"
        if os.path.exists(os.path.join(output_dir, montage_path)):
            html_content += f"""
            <div class="plot-container">
                <h3>ROI: {roi}</h3>
                <img src="{montage_path}" alt="Montage Plot for {roi}">
            </div>
            """
            
    html_content += """
    </body>
    </html>
    """
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
        
    logger.info(f"HTML report generation complete: {report_path}")

def main(config_path='config/config.yaml'):
    logger.info("========== Step 3b: Visualization and Reporting ==========")
    config = parsing.load_yaml(config_path)
    if config is None: return 

    try:
        output_dir = config['output_dir']
        rois = [r['name'] for r in config['roi_info']]
        groups = [g['name'] for g in config['group_info']]
    except KeyError as e:
        logger.error(f"config.yaml key error: {e}")
        return
        
    agg_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    main_stats_path = os.path.join(output_dir, "statistical_results_main.csv")
    posthoc_stats_path = os.path.join(output_dir, "statistical_results_posthoc.csv")
    
    if not os.path.exists(agg_csv_path):
        logger.error("Aggregated file not found. Please run Step 2 first.")
        return

    # Load data
    df_agg = pd.read_csv(agg_csv_path)
    id_vars = ['group', 'n', 'roi']
    value_vars = [col for col in df_agg.columns if col not in id_vars]
    df_long = df_agg.melt(id_vars=id_vars, value_vars=value_vars, var_name='m_z_bin', value_name='intensity')
    
    try:
        df_long['group'] = pd.Categorical(df_long['group'], categories=groups, ordered=True)
    except Exception as e:
        logger.warning(f"Group order application warning: {e}")

    df_main_stats = pd.DataFrame()
    if os.path.exists(main_stats_path):
        df_main_stats = pd.read_csv(main_stats_path)
        
    df_posthoc_stats = pd.DataFrame()
    if os.path.exists(posthoc_stats_path):
        df_posthoc_stats = pd.read_csv(posthoc_stats_path)

    palette = sns.color_palette("Paired", n_colors=len(groups))
    group_color_map = dict(zip(groups, palette))
    
    for roi in rois:
        logger.info(f"Visualizing ROI: {roi}")
        df_long_roi = df_long[df_long['roi'] == roi].copy()
        if df_long_roi.empty:
            logger.warning(f"No data: {roi}")
            continue
            
        df_posthoc_roi = pd.DataFrame()
        if not df_posthoc_stats.empty:
             df_posthoc_roi = df_posthoc_stats[df_posthoc_stats['roi'] == roi]

        # Generate individual plots
        for m_z_bin in value_vars:
            data_bin = df_long_roi[df_long_roi['m_z_bin'] == m_z_bin]
            stats_df_bin = pd.DataFrame()
            if not df_posthoc_roi.empty:
                stats_df_bin = df_posthoc_roi[df_posthoc_roi['m_z_bin'] == m_z_bin]
            generate_single_plot(data_bin, roi, m_z_bin, stats_df_bin, group_color_map, groups, output_dir)

        # Generate montage plot
        generate_montage_plot(df_long_roi, roi, value_vars, df_posthoc_roi, group_color_map, groups, output_dir)

    # Generate HTML Report
    generate_html_report(config, output_dir, df_main_stats, df_posthoc_stats, rois)

    logger.info("========== Step 3b complete ==========")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Step 3b: Visualization and Reporting',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config/config.yaml',
        help='Path to config YAML file'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help='Path to log file'
    )
    
    args = parser.parse_args()
    
    global logger
    logger = get_logger(__name__, log_file=args.log_file)
    
    main(config_path=args.config)
