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
    # Robust check for significance (handles boolean True, string 'True', 'true', etc.)
    if not stats_df_bin.empty and 'significant' in stats_df_bin.columns:
        # Ensure boolean type for filtering
        is_sig = stats_df_bin['significant'].astype(str).str.lower() == 'true'
        significant_pairs = stats_df_bin[is_sig]
    else:
        significant_pairs = pd.DataFrame()
    
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
    # Limit columns to max 4 to prevent plots from becoming too small
    cols = min(4, int(np.ceil(np.sqrt(num_bins))))
    if cols < 1: cols = 1 # Safety check
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
    logger.info("Starting HTML report generation (Enhanced)")
    
    report_path = os.path.join(output_dir, 'analysis_report.html')
    
    # Calculate summary stats
    n_rois = len(rois)
    n_significant = 0
    if not df_main_stats.empty and 'significant' in df_main_stats.columns:
        n_significant = df_main_stats['significant'].astype(str).str.lower().eq('true').sum()
        
    generated_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # Extract config info for dashboard
    groups = [g['name'] for g in config.get('group_info', [])]
    group_str = ", ".join(groups)
    
    # Get m/z bins info (if available from stats or config)
    mz_bins = []
    if not df_main_stats.empty and 'm_z_bin' in df_main_stats.columns:
        mz_bins = df_main_stats['m_z_bin'].unique().tolist()
    n_bins = len(mz_bins)
    
    # Extract settings for Config tab
    outlier_settings = config.get('outlier_detection', {})
    stats_settings = config.get('statistics_settings', {})
    binning_settings = config.get('binning_settings', {})

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Mass Imaging Analysis Report</title>
        
        <!-- Bootstrap CSS -->
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <!-- DataTables CSS -->
        <link href="https://cdn.datatables.net/1.13.4/css/dataTables.bootstrap5.min.css" rel="stylesheet">
        
        <style>
            body {{ background-color: #f8f9fa; padding-top: 20px; }}
            .container {{ max-width: 1400px; }}
            .card {{ margin-bottom: 20px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border: none; }}
            .card-header {{ background-color: #fff; border-bottom: 1px solid #eee; font-weight: bold; }}
            .plot-img {{ width: 100%; height: auto; border-radius: 4px; transition: transform 0.2s; }}
            .plot-img:hover {{ transform: scale(1.02); }}
            .nav-tabs .nav-link {{ color: #495057; }}
            .nav-tabs .nav-link.active {{ font-weight: bold; color: #0d6efd; }}
            .summary-box {{ padding: 20px; background: white; border-radius: 8px; height: 100%; }}
            .summary-title {{ font-size: 0.9rem; color: #6c757d; text-transform: uppercase; letter-spacing: 1px; margin-bottom: 10px; }}
            .summary-content {{ font-size: 1.1rem; font-weight: 500; color: #212529; }}
            .summary-highlight {{ color: #0d6efd; font-weight: bold; }}
            pre {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #e9ecef; }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="d-flex justify-content-between align-items-center mb-4">
                <div>
                    <h1 class="display-5 fw-bold">Mass Imaging Analysis Report</h1>
                    <p class="text-muted">Generated on {generated_time}</p>
                </div>
                <div>
                    <span class="badge bg-primary rounded-pill">v1.1</span>
                </div>
            </div>

            <!-- Summary Dashboard -->
            <div class="row mb-4">
                <div class="col-md-3">
                    <div class="summary-box">
                        <div class="summary-title">Groups ({len(groups)})</div>
                        <div class="summary-content">{group_str}</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="summary-box">
                        <div class="summary-title">ROIs ({n_rois})</div>
                        <div class="summary-content">{", ".join(rois)}</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="summary-box">
                        <div class="summary-title">m/z Bins ({n_bins})</div>
                        <div class="summary-content">Analyzed {n_bins} mass bins</div>
                    </div>
                </div>
                <div class="col-md-3">
                    <div class="summary-box">
                        <div class="summary-title">Significance</div>
                        <div class="summary-content"><span class="summary-highlight">{n_significant}</span> significant tests found</div>
                    </div>
                </div>
            </div>

            <!-- Tabs -->
            <ul class="nav nav-tabs mb-4" id="reportTabs" role="tablist">
                <li class="nav-item" role="presentation">
                    <button class="nav-link active" id="plots-tab" data-bs-toggle="tab" data-bs-target="#plots" type="button" role="tab">Plots</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="tables-tab" data-bs-toggle="tab" data-bs-target="#tables" type="button" role="tab">Statistical Tables</button>
                </li>
                <li class="nav-item" role="presentation">
                    <button class="nav-link" id="config-tab" data-bs-toggle="tab" data-bs-target="#config" type="button" role="tab">Configuration</button>
                </li>
            </ul>

            <div class="tab-content" id="reportTabsContent">
                
                <!-- Plots Tab -->
                <div class="tab-pane fade show active" id="plots" role="tabpanel">
                    <div class="row">
    """
    
    for roi in rois:
        montage_path = f"plot_montage_roi_{roi}.png"
        if os.path.exists(os.path.join(output_dir, montage_path)):
            html_content += f"""
                        <div class="col-md-12 mb-4">
                            <div class="card">
                                <div class="card-header d-flex justify-content-between align-items-center">
                                    <span>ROI: {roi}</span>
                                    <a href="{montage_path}" target="_blank" class="btn btn-sm btn-outline-primary">Open Full Size</a>
                                </div>
                                <div class="card-body">
                                    <img src="{montage_path}" class="plot-img" alt="Montage Plot for {roi}">
                                </div>
                            </div>
                        </div>
            """
            
    html_content += """
                    </div>
                </div>

                <!-- Tables Tab -->
                <div class="tab-pane fade" id="tables" role="tabpanel">
                    <div class="card mb-4">
                        <div class="card-header">Main Statistical Results</div>
                        <div class="card-body">
                            <div class="table-responsive">
                                <table id="mainTable" class="table table-striped table-hover" style="width:100%">
                                    <thead>
                                        <tr>
    """
    
    if not df_main_stats.empty:
        for col in df_main_stats.columns:
            html_content += f"<th>{col}</th>"
        html_content += "</tr></thead><tbody>"
        
        for _, row in df_main_stats.iterrows():
            html_content += "<tr>"
            for val in row:
                # Format floats
                if isinstance(val, (float, np.floating)):
                    html_content += f"<td>{val:.4e}</td>"
                else:
                    html_content += f"<td>{val}</td>"
            html_content += "</tr>"
        html_content += "</tbody></table></div></div></div>"
    else:
        html_content += "<th>No Data</th></tr></thead><tbody></tbody></table></div></div></div>"

    if not df_posthoc_stats.empty:
        html_content += """
                    <div class="card mb-4">
                        <div class="card-header">Post-hoc Results</div>
                        <div class="card-body">
                            <div class="table-responsive">
                                <table id="posthocTable" class="table table-striped table-hover" style="width:100%">
                                    <thead>
                                        <tr>
        """
        for col in df_posthoc_stats.columns:
            html_content += f"<th>{col}</th>"
        html_content += "</tr></thead><tbody>"
        
        for _, row in df_posthoc_stats.iterrows():
            html_content += "<tr>"
            for val in row:
                if isinstance(val, (float, np.floating)):
                    html_content += f"<td>{val:.4e}</td>"
                else:
                    html_content += f"<td>{val}</td>"
            html_content += "</tr>"
        html_content += "</tbody></table></div></div></div>"

    # Helper to format dict as string
    def format_dict(d):
        return "\n".join([f"{k}: {v}" for k, v in d.items()])

    html_content += f"""
                </div>

                <!-- Config Tab -->
                <div class="tab-pane fade" id="config" role="tabpanel">
                    <div class="row">
                        <div class="col-md-6">
                            <div class="card">
                                <div class="card-header">General Settings</div>
                                <div class="card-body">
                                    <strong>Output Directory:</strong> {output_dir}<br>
                                    <strong>Binning Mode:</strong> {binning_settings.get('mode', 'N/A')}<br>
                                    <strong>Binning File:</strong> {binning_settings.get('file_path', 'N/A')}
                                </div>
                            </div>
                        </div>
                        <div class="col-md-6">
                            <div class="card">
                                <div class="card-header">Statistics Settings</div>
                                <div class="card-body">
                                    <pre>{format_dict(stats_settings)}</pre>
                                </div>
                            </div>
                        </div>
                        <div class="col-md-6">
                            <div class="card">
                                <div class="card-header">Outlier Detection</div>
                                <div class="card-body">
                                    <pre>{format_dict(outlier_settings)}</pre>
                                </div>
                            </div>
                        </div>
                        <div class="col-md-6">
                            <div class="card">
                                <div class="card-header">Groups & ROIs</div>
                                <div class="card-body">
                                    <strong>Groups:</strong> {group_str}<br>
                                    <strong>ROIs:</strong> {", ".join(rois)}
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Scripts -->
        <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.4/js/dataTables.bootstrap5.min.js"></script>
        <script>
            $(document).ready(function() {{
                $('#mainTable').DataTable({{
                    "pageLength": 10,
                    "order": [[ 3, "asc" ]] 
                }});
                $('#posthocTable').DataTable({{
                    "pageLength": 10
                }});
            }});
        </script>
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
        # Ensure m_z_bin is string for consistent matching
        if 'm_z_bin' in df_posthoc_stats.columns:
            df_posthoc_stats['m_z_bin'] = df_posthoc_stats['m_z_bin'].astype(str)

    # Ensure m_z_bin in df_long is also string
    df_long['m_z_bin'] = df_long['m_z_bin'].astype(str)
    value_vars = [str(v) for v in value_vars] # Update value_vars to strings too

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
