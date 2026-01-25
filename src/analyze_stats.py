import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')
warnings.filterwarnings('ignore', category=FutureWarning, module='seaborn')

import os
import pandas as pd
import itertools
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

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

def run_statistics(data_long_roi, m_z_bin, roi, test_type, p_threshold):
    data = data_long_roi[data_long_roi['m_z_bin'] == m_z_bin]
    groups = data['group'].unique()
    n_groups = len(groups)
    
    main_test_result = {
        'roi': roi, 
        'm_z_bin': m_z_bin, 
        'test_name': 'N/A', 
        'p_value': 1.0, 
        'significant': False
    }
    posthoc_results = [] 

    if n_groups < 2:
        return main_test_result, posthoc_results

    all_group_data = [data[data['group'] == g]['intensity'] for g in groups]

    try:
        if n_groups == 2:
            g1_data = all_group_data[0]
            g2_data = all_group_data[1]
            g1_name = groups[0]
            g2_name = groups[1]

            if test_type == 'parametric':
                test_name = 't-test_ind'
                stat, p_val = ttest_ind(g1_data, g2_data, equal_var=False) 
            else:
                test_name = 'Mann-Whitney U'
                if g1_data.empty or g2_data.empty or np.array_equal(g1_data.values, g2_data.values):
                     p_val = 1.0
                else:
                    stat, p_val = mannwhitneyu(g1_data, g2_data, alternative='two-sided')
            
            main_test_result['test_name'] = test_name
            main_test_result['p_value'] = p_val
            main_test_result['significant'] = p_val < p_threshold
            
            posthoc_results.append({
                'roi': roi,
                'm_z_bin': m_z_bin,
                'test_name': test_name,
                'group1': g1_name,
                'group2': g2_name,
                'p_value': p_val,
                'p_adj': p_val, 
                'significant': p_val < p_threshold
            })

        else:
            if test_type == 'parametric':
                test_name = 'ANOVA'
                stat, p_val = f_oneway(*all_group_data)
            else:
                test_name = 'Kruskal-Wallis'
                stat, p_val = kruskal(*all_group_data)
            
            main_test_result['test_name'] = test_name
            main_test_result['p_value'] = p_val
            main_test_result['significant'] = p_val < p_threshold

            if main_test_result['significant']:
                if test_type == 'parametric':
                    # Tukey's HSD
                    df_tukey = data[['intensity', 'group']]
                    tukey_result = pairwise_tukeyhsd(df_tukey['intensity'], df_tukey['group'], alpha=p_threshold)
                    df_posthoc = pd.DataFrame(data=tukey_result._results_table.data[1:], columns=tukey_result._results_table.data[0])
                    for _, row in df_posthoc.iterrows():
                        posthoc_results.append({
                            'roi': roi, 'm_z_bin': m_z_bin, 'test_name': 'Tukey HSD',
                            'group1': row['group1'], 'group2': row['group2'],
                            'p_value': row['p-adj'], 'p_adj': row['p-adj'], 'significant': row['reject'] 
                        })
                else:
                    # Bonferroni Mann-Whitney
                    posthoc_test_name = 'Mann-Whitney (Bonferroni)'
                    all_pairs = list(itertools.combinations(groups, 2))
                    num_comparisons = len(all_pairs)
                    for g1_name, g2_name in all_pairs:
                        g1_data = data[data['group'] == g1_name]['intensity']
                        g2_data = data[data['group'] == g2_name]['intensity']
                        if g1_data.empty or g2_data.empty or np.array_equal(g1_data.values, g2_data.values):
                            p_val = 1.0
                            p_adj = 1.0
                        else:
                            stat, p_val = mannwhitneyu(g1_data, g2_data, alternative='two-sided')
                            p_adj = min(p_val * num_comparisons, 1.0)
                        posthoc_results.append({
                            'roi': roi, 'm_z_bin': m_z_bin, 'test_name': posthoc_test_name,
                            'group1': g1_name, 'group2': g2_name,
                            'p_value': p_val, 'p_adj': p_adj, 'significant': p_adj < p_threshold
                        })
                        
    except Exception as e:
        logger.error(f"Statistical analysis error (bin: {m_z_bin}): {e}", exc_info=True)

    return main_test_result, posthoc_results

def apply_correction(stats_results, method, p_threshold):
    """
    Apply multiple comparison correction to main statistical results.
    
    Args:
        stats_results: List of main test result dictionaries
        method: 'fdr_bh', 'bonferroni', or 'none'
        p_threshold: Original p-value threshold
    
    Returns:
        List of corrected result dictionaries with 'p_adj' and updated 'significant'
    """
    if not stats_results or len(stats_results) == 0:
        return stats_results
    
    if method == 'none':
        # No correction - add p_adj = p_value
        for result in stats_results:
            result['p_adj'] = result['p_value']
        return stats_results
    
    # Extract p-values
    p_values = [r['p_value'] for r in stats_results]
    
    try:
        if method == 'fdr_bh':
            # Benjamini-Hochberg FDR
            reject, p_adj, _, _ = multipletests(p_values, alpha=p_threshold, method='fdr_bh')
        elif method == 'bonferroni':
            # Bonferroni correction
            reject, p_adj, _, _ = multipletests(p_values, alpha=p_threshold, method='bonferroni')
        else:
            logger.warning(f"Unknown correction method '{method}', using 'none'")
            for result in stats_results:
                result['p_adj'] = result['p_value']
            return stats_results
        
        # Update results with adjusted p-values
        for i, result in enumerate(stats_results):
            result['p_adj'] = p_adj[i]
            result['significant'] = reject[i]
        
    except Exception as e:
        logger.error(f"Error applying correction: {e}")
        # Fallback: no correction
        for result in stats_results:
            result['p_adj'] = result['p_value']
    
    return stats_results

def get_control_group(config):
    """
    Get the control/reference group name from config.
    
    Args:
        config: Configuration dictionary
    
    Returns:
        str: Name of control group
    """
    groups_info = config['group_info']
    
    # Find group with is_control: true
    control_groups = [g['name'] for g in groups_info if g.get('is_control', False)]
    
    if len(control_groups) == 0:
        # Default to first group if no control specified
        logger.warning("No control group specified (is_control: true), using first group as control")
        return groups_info[0]['name']
    elif len(control_groups) > 1:
        logger.warning(f"Multiple control groups specified: {control_groups}, using first one")
        return control_groups[0]
    else:
        return control_groups[0]

def create_integrated_results(df_agg, df_stats, config, roi):
    """
    Create DESeq2-style integrated results table combining:
    - Individual sample intensities
    - Group statistics (mean, SD)
    - Fold changes (log2 and raw)
    - Statistical test results
    
    Args:
        df_agg: Aggregated mean intensities DataFrame (wide format)
        df_stats: Statistical results DataFrame
        config: Configuration dictionary
        roi: ROI name
    
    Returns:
        DataFrame with integrated results
    """
    try:
        # Get control group
        control_group = get_control_group(config)
        groups = [g['name'] for g in config['group_info']]
        
        logger.info(f"Creating integrated results for ROI: {roi}, control group: {control_group}")
        
        # Filter for this ROI
        df_agg_roi = df_agg[df_agg['roi'] == roi].copy()
        df_stats_roi = df_stats[df_stats['roi'] == roi].copy()
        
        if df_agg_roi.empty:
            logger.warning(f"No aggregated data for ROI: {roi}")
            return pd.DataFrame()
        
        # Get m/z bins
        id_vars = ['group', 'n', 'roi']
        m_z_bins = [col for col in df_agg_roi.columns if col not in id_vars]
        
        # Convert to long format for calculations
        df_long = df_agg_roi.melt(
            id_vars=id_vars,
            value_vars=m_z_bins,
            var_name='m_z_bin',
            value_name='intensity'
        )
        
        # Calculate group statistics
        group_stats_list = []
        for m_z_bin in m_z_bins:
            df_bin = df_long[df_long['m_z_bin'] == m_z_bin]
            
            # Overall mean
            base_mean = df_bin['intensity'].mean()
            
            # Per-group statistics
            stats_row = {
                'roi': roi,
                'm_z_bin': m_z_bin,
                'baseMean': base_mean
            }
            
            group_means = {}
            for group in groups:
                df_group = df_bin[df_bin['group'] == group]
                group_mean = df_group['intensity'].mean()
                group_sd = df_group['intensity'].std()
                
                stats_row[f'{group}_mean'] = group_mean
                stats_row[f'{group}_sd'] = group_sd
                group_means[group] = group_mean
            
            # Calculate fold changes vs control
            control_mean = group_means.get(control_group, 1.0)
            if control_mean == 0:
                control_mean = 1e-10  # Avoid division by zero
            
            # Calculate fold changes for each non-control group
            for group in groups:
                if group != control_group:
                    experimental_mean = group_means[group]
                    fold_change = experimental_mean / control_mean
                    log2_fc = np.log2(fold_change) if fold_change > 0 else np.nan
                    
                    stats_row[f'log2FC_{group}_vs_{control_group}'] = log2_fc
                    stats_row[f'FC_{group}_vs_{control_group}'] = fold_change
            
            group_stats_list.append(stats_row)
        
        df_group_stats = pd.DataFrame(group_stats_list)
        
        # Merge with statistical results
        df_integrated = df_group_stats.merge(
            df_stats_roi[['m_z_bin', 'test_name', 'p_value', 'p_adj', 'significant']],
            on='m_z_bin',
            how='left'
        )
        
        # Add individual sample data
        for group in groups:
            df_group = df_agg_roi[df_agg_roi['group'] == group].copy()
            n_per_group = df_group['n'].max()
            
            for n in range(1, int(n_per_group) + 1):
                df_sample = df_group[df_group['n'] == n]
                if not df_sample.empty:
                    # Get sample values for each m/z bin
                    sample_col_name = f'{group}_{n}'
                    sample_values = df_sample[m_z_bins].iloc[0].to_dict()
                    
                    # Add to integrated dataframe
                    for m_z_bin in m_z_bins:
                        mask = df_integrated['m_z_bin'] == m_z_bin
                        df_integrated.loc[mask, sample_col_name] = sample_values.get(m_z_bin, np.nan)
        
        # Reorder columns for better readability
        # Order: identifiers, summary stats, fold changes, statistics, samples
        base_cols = ['roi', 'm_z_bin', 'baseMean']
        group_mean_cols = [f'{g}_mean' for g in groups]
        group_sd_cols = [f'{g}_sd' for g in groups]
        
        fc_cols = []
        for group in groups:
            if group != control_group:
                fc_cols.append(f'log2FC_{group}_vs_{control_group}')
                fc_cols.append(f'FC_{group}_vs_{control_group}')
        
        stat_cols = ['pvalue', 'p_adj', 'significant', 'test_name']
        # Rename for consistency
        df_integrated = df_integrated.rename(columns={'p_value': 'pvalue'})
        
        sample_cols = []
        for group in groups:
            n_per_group = len(df_agg_roi[df_agg_roi['group'] == group])
            for n in range(1, n_per_group + 1):
                col_name = f'{group}_{n}'
                if col_name in df_integrated.columns:
                    sample_cols.append(col_name)
        
        # Combine in desired order
        ordered_cols = base_cols + group_mean_cols + group_sd_cols + fc_cols + stat_cols + sample_cols
        
        # Only include columns that exist
        ordered_cols = [col for col in ordered_cols if col in df_integrated.columns]
        df_integrated = df_integrated[ordered_cols]
        
        logger.info(f"Created integrated results: {len(df_integrated)} rows, {len(df_integrated.columns)} columns")
        
        return df_integrated
        
    except Exception as e:
        logger.error(f"Error creating integrated results for ROI {roi}: {e}", exc_info=True)
        return pd.DataFrame()

def export_to_prism(df_long_roi, roi, config, m_z_bins, output_dir):
    try:
        logger.info(f"Generating Prism CSV: {roi}")
        groups_info_list = config['group_info']
        prism_columns = []
        for group_dict in groups_info_list:
            g = group_dict['name']
            for n in range(1, group_dict['n_per_group'] + 1):
                prism_columns.append(f"{g}-{n}")
        
        prism_data = {} 
        for m_z_bin in m_z_bins:
            row_data = []
            for group_dict in groups_info_list:
                g = group_dict['name']
                for n in range(1, group_dict['n_per_group'] + 1):
                    val_series = df_long_roi[(df_long_roi['m_z_bin'] == m_z_bin) & (df_long_roi['group'] == g) & (df_long_roi['n'] == n)]['intensity']
                    row_data.append(val_series.iloc[0] if not val_series.empty else None)
            prism_data[m_z_bin] = row_data

        df_prism = pd.DataFrame.from_dict(prism_data, orient='index', columns=prism_columns)
        prism_file_path = os.path.join(output_dir, f"aggregated_data_roi_{roi}_prism.csv")
        df_prism.to_csv(prism_file_path, float_format='%.4f', index=True, index_label='m_z_bin')
        
    except Exception as e:
        logger.error(f"Prism conversion error: {e}", exc_info=True)

def main(config_path='config/config.yaml'):
    logger.info("========== Step 3a: Starting statistical analysis ==========")
    config = parsing.load_yaml(config_path)
    if config is None: return 

    try:
        output_dir = config['output_dir']
        rois = [r['name'] for r in config['roi_info']]
        stats_settings = config['statistics_settings']
        test_type = stats_settings['test_type']
        p_threshold = stats_settings['p_value_threshold']
        correction_method = stats_settings.get('multiple_comparison_correction', 'fdr_bh')
        groups = [g['name'] for g in config['group_info']]
        
        logger.info(f"Multiple comparison correction: {correction_method}")
    except KeyError as e:
        logger.error(f"config.yaml key error: {e}")
        return
        
    agg_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    if not os.path.exists(agg_csv_path):
        logger.error("Aggregated file not found. Please run Step 2 first.")
        return
        
    df_agg = pd.read_csv(agg_csv_path)
    id_vars = ['group', 'n', 'roi']
    value_vars = [col for col in df_agg.columns if col not in id_vars]
    
    df_long = df_agg.melt(id_vars=id_vars, value_vars=value_vars, var_name='m_z_bin', value_name='intensity')
    
    try:
        df_long['group'] = pd.Categorical(df_long['group'], categories=groups, ordered=True)
    except Exception as e:
        logger.warning(f"Group order application warning: {e}")
    
    all_main_stats = []
    all_posthoc_stats = []
    
    for roi in rois:
        logger.info(f"Analyzing ROI: {roi}")
        df_long_roi = df_long[df_long['roi'] == roi].copy()
        if df_long_roi.empty:
            logger.warning(f"No data: {roi}")
            continue
        
        # Collect results for this ROI
        roi_main_stats = []
        roi_posthoc_stats = []
        
        for m_z_bin in value_vars: 
            main_result, posthoc_results = run_statistics(df_long_roi, m_z_bin, roi, test_type, p_threshold)
            roi_main_stats.append(main_result)
            roi_posthoc_stats.extend(posthoc_results)
        
        # Apply multiple comparison correction to main tests for this ROI
        roi_main_stats = apply_correction(roi_main_stats, correction_method, p_threshold)
        logger.info(f"Applied {correction_method} correction to {len(roi_main_stats)} tests for ROI: {roi}")
        
        # Add to overall results
        all_main_stats.extend(roi_main_stats)
        all_posthoc_stats.extend(roi_posthoc_stats)
        
        # Export to Prism format
        export_to_prism(df_long_roi, roi, config, value_vars, output_dir)
    
    # Create integrated results tables
    logger.info("========== Creating integrated results tables ==========")
    all_integrated_results = []
    
    # Convert main stats list to DataFrame for merging
    df_main_stats = pd.DataFrame(all_main_stats)
    
    for roi in rois:
        logger.info(f"Creating integrated results for ROI: {roi}")
        
        # Create integrated results for this ROI
        df_integrated_roi = create_integrated_results(
            df_agg, df_main_stats, config, roi
        )
        
        if not df_integrated_roi.empty:
            # Save per-ROI integrated results
            integrated_roi_path = os.path.join(output_dir, f"integrated_results_{roi}.csv")
            df_integrated_roi.to_csv(integrated_roi_path, index=False, float_format='%.4f')
            logger.info(f"Saved integrated results for {roi}: {integrated_roi_path}")
            
            all_integrated_results.append(df_integrated_roi)
    
    # Save combined integrated results
    if all_integrated_results:
        df_integrated_all = pd.concat(all_integrated_results, ignore_index=True)
        integrated_all_path = os.path.join(output_dir, "integrated_results_all.csv")
        df_integrated_all.to_csv(integrated_all_path, index=False, float_format='%.4f')
        logger.info(f"Saved combined integrated results: {integrated_all_path}")
        logger.info(f"Total integrated results: {len(df_integrated_all)} rows across {len(rois)} ROIs")

    try:
        df_main_stats = pd.DataFrame(all_main_stats)
        main_stats_path = os.path.join(output_dir, "statistical_results_main.csv")
        df_main_stats.to_csv(main_stats_path, index=False, float_format='%.4e')
        
        df_posthoc_stats = pd.DataFrame()
        if all_posthoc_stats: 
            df_posthoc_stats = pd.DataFrame(all_posthoc_stats)
            posthoc_stats_path = os.path.join(output_dir, "statistical_results_posthoc.csv")
            df_posthoc_stats.to_csv(posthoc_stats_path, index=False, float_format='%.4e')
        else:
            logger.info("No post-hoc results (no significance).")

    except Exception as e:
        logger.error(f"Error saving results: {e}", exc_info=True)

    logger.info("========== Step 3a complete ==========")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Step 3a: Perform statistical analysis',
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
