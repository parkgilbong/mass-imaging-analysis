import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')
warnings.filterwarnings('ignore', category=FutureWarning, module='seaborn')

import os
import pandas as pd
import itertools
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd

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
        groups = [g['name'] for g in config['group_info']]
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
        
        for m_z_bin in value_vars: 
            main_result, posthoc_results = run_statistics(df_long_roi, m_z_bin, roi, test_type, p_threshold)
            all_main_stats.append(main_result)
            all_posthoc_stats.extend(posthoc_results) 
        
        # Export to Prism format
        export_to_prism(df_long_roi, roi, config, value_vars, output_dir)

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
