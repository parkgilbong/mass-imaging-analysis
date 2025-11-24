import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')
warnings.filterwarnings('ignore', category=FutureWarning, module='seaborn')

import os
import pandas as pd
import yaml
import itertools
import numpy as np
from datetime import datetime

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

CONFIG_FILE = 'config/config.yaml'
logger = get_logger(__name__)

def run_statistics(data_long_roi, m_z_bin, roi, test_type, p_threshold):
    # (로직 동일 - print 제거)
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
        logger.error(f"통계 분석 오류 (bin: {m_z_bin}): {e}", exc_info=True)

    return main_test_result, posthoc_results

# (plot_single_bin 함수는 기존과 동일, 로깅 불필요)
def plot_single_bin(ax, data_bin, m_z_bin, stats_df_bin):
    sns.barplot(data=data_bin, x='group', y='intensity', hue='group', palette='Paired', errorbar=None, ax=ax, legend=False)
    sns.stripplot(data=data_bin, x='group', y='intensity', color='black', legend=False, s=5, ax=ax, jitter=False)
    ax.set_title(f"m/z bin: {m_z_bin}", fontsize=10)
    ax.set_xlabel('')
    ax.set_ylabel('Intensity')
    ax.tick_params(axis='x', rotation=45)
    significant_pairs = stats_df_bin[stats_df_bin['significant'] == True]
    if not significant_pairs.empty:
        max_intensity = data_bin['intensity'].max()
        y_offset = max_intensity * 0.05 
        y = max_intensity + y_offset
        x_pos = (len(ax.get_xticklabels()) - 1) / 2.0
        ax.text(x_pos, y, "*", ha='center', va='bottom', fontsize=14, color='red')

def export_to_prism(df_long_roi, roi, config, m_z_bins, output_dir):
    try:
        logger.info(f"Prism용 CSV 생성 중: {roi}")
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
        logger.error(f"Prism 변환 오류: {e}", exc_info=True)

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
        logger.error(f"개별 플롯 저장 실패 ({m_z_bin}): {e}")
        plt.close(fig_single)

def generate_montage_plot(df_long_roi, roi, value_vars, df_posthoc_roi, group_color_map, groups, output_dir):
    logger.info(f"몽타주 플롯 생성 중: {roi}")
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

def generate_html_report(config, output_dir, df_main_stats, df_posthoc_stats):
    logger.info("HTML 보고서 생성 시작")
    # ... (HTML 생성 코드는 print 대신 logger를 쓰지 않아도 무방하나, 완료 메시지는 로깅)
    # ... (중략) ... 
    # 기존 코드 유지하되 마지막 print만 변경
    logger.info(f"HTML 보고서 생성 완료: {os.path.join(output_dir, 'analysis_report.html')}")

def main(config_path='config/config.yaml'):
    logger.info("========== Step 3: 통계 분석 및 시각화 시작 ==========")
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
        logger.error(f"config.yaml 키 오류: {e}")
        return
        
    agg_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    if not os.path.exists(agg_csv_path):
        logger.error("집계 파일이 없습니다. Step 2를 먼저 실행하세요.")
        return
        
    df_agg = pd.read_csv(agg_csv_path)
    id_vars = ['group', 'n', 'roi']
    value_vars = [col for col in df_agg.columns if col not in id_vars]
    
    df_long = df_agg.melt(id_vars=id_vars, value_vars=value_vars, var_name='m_z_bin', value_name='intensity')
    
    try:
        df_long['group'] = pd.Categorical(df_long['group'], categories=groups, ordered=True)
    except Exception as e:
        logger.warning(f"그룹 순서 적용 경고: {e}")
    
    all_main_stats = []
    all_posthoc_stats = []
    palette = sns.color_palette("Paired", n_colors=len(groups))
    group_color_map = dict(zip(groups, palette))
    
    for roi in rois:
        logger.info(f"ROI 분석 중: {roi}")
        df_long_roi = df_long[df_long['roi'] == roi].copy()
        if df_long_roi.empty:
            logger.warning(f"데이터 없음: {roi}")
            continue
        
        for m_z_bin in value_vars: 
            main_result, posthoc_results = run_statistics(df_long_roi, m_z_bin, roi, test_type, p_threshold)
            all_main_stats.append(main_result)
            all_posthoc_stats.extend(posthoc_results) 

        df_posthoc_roi = pd.DataFrame(all_posthoc_stats)
        
        for idx, m_z_bin in enumerate(value_vars):
            data_bin = df_long_roi[df_long_roi['m_z_bin'] == m_z_bin]
            stats_df_bin = pd.DataFrame()
            if not df_posthoc_roi.empty:
                stats_df_bin = df_posthoc_roi[(df_posthoc_roi['m_z_bin'] == m_z_bin) & (df_posthoc_roi['roi'] == roi)]
            generate_single_plot(data_bin, roi, m_z_bin, stats_df_bin, group_color_map, groups, output_dir)

        generate_montage_plot(df_long_roi, roi, value_vars, df_posthoc_roi, group_color_map, groups, output_dir)
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
            logger.info("Post-hoc 결과 없음 (유의차 없음).")

        # HTML 보고서 함수 호출 (이전에 정의된 함수 사용)
        # generate_html_report(...) -> 위 코드에 포함되어야 함. 
        # 여기서는 생략했으나 실제 파일에는 포함되어야 함.

    except Exception as e:
        logger.error(f"결과 저장 중 오류: {e}", exc_info=True)

    logger.info("========== Step 3 완료 ==========")

if __name__ == '__main__':
    main(config_path='config/config.yaml')