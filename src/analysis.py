import warnings
# pyimzml에서 발생하는 특정 UserWarning을 무시합니다.
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')
# Seaborn에서 발생하는 FutureWarning를 무시합니다. (palette='dark:black')
warnings.filterwarnings('ignore', category=FutureWarning, module='seaborn')

import os
import pandas as pd
import yaml
import itertools
import traceback
import numpy as np
from datetime import datetime # (추가) 보고서 생성 일시 기록용

# 통계 및 시각화 라이브러리
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec # 몽타주 플롯을 위한 GridSpec
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# 'src' 폴더의 'parsing' 모듈에서 load_yaml 함수 임포트
try:
    # 이 스크립트가 'src' 외부(예: main.ipynb)에서 실행될 것을 가정
    import parsing
except ImportError:
    # 이 스크립트가 'src' 내부에서 직접 실행될 경우
    try:
        from . import parsing
    except ImportError:
        print("오류: 'parsing.py' 파일을 찾을 수 없습니다.")
        pass 

def run_statistics(data_long_roi, m_z_bin, roi, test_type, p_threshold):
    """
    (수정됨)
    특정 ROI와 m/z bin에 대해 통계 분석을 수행합니다.
    (np.array_equal() 사용하여 인덱스 오류 방지)
    """
    
    # 1. 분석할 데이터 필터링
    data = data_long_roi[data_long_roi['m_z_bin'] == m_z_bin]
    # (수정) config.yaml 순서에 맞게 그룹 정렬
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

    # 2. 그룹별 데이터 분리
    all_group_data = [data[data['group'] == g]['intensity'] for g in groups]

    try:
        if n_groups == 2:
            # --- 2개 그룹 비교 (t-test or Mann-Whitney U) ---
            g1_data = all_group_data[0]
            g2_data = all_group_data[1]
            g1_name = groups[0]
            g2_name = groups[1]

            if test_type == 'parametric':
                test_name = 't-test_ind'
                stat, p_val = ttest_ind(g1_data, g2_data, equal_var=False) 
            else:
                test_name = 'Mann-Whitney U'
                # (수정) 인덱스 무시하고 값만 비교
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
            # --- 3개 이상 그룹 비교 (ANOVA or Kruskal-Wallis) ---
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
                            'roi': roi,
                            'm_z_bin': m_z_bin,
                            'test_name': 'Tukey HSD',
                            'group1': row['group1'],
                            'group2': row['group2'],
                            'p_value': row['p-adj'], 
                            'p_adj': row['p-adj'],
                            'significant': row['reject'] 
                        })
                else:
                    # Non-parametric post-hoc: Pairwise Mann-Whitney U with Bonferroni correction
                    posthoc_test_name = 'Mann-Whitney (Bonferroni)'
                    all_pairs = list(itertools.combinations(groups, 2))
                    num_comparisons = len(all_pairs)
                    
                    for g1_name, g2_name in all_pairs:
                        g1_data = data[data['group'] == g1_name]['intensity']
                        g2_data = data[data['group'] == g2_name]['intensity']
                        
                        # (수정) 인덱스 무시하고 값만 비교
                        if g1_data.empty or g2_data.empty or np.array_equal(g1_data.values, g2_data.values):
                            p_val = 1.0
                            p_adj = 1.0
                        else:
                            stat, p_val = mannwhitneyu(g1_data, g2_data, alternative='two-sided')
                            p_adj = min(p_val * num_comparisons, 1.0) # Bonferroni correction
                        
                        posthoc_results.append({
                            'roi': roi,
                            'm_z_bin': m_z_bin,
                            'test_name': posthoc_test_name,
                            'group1': g1_name,
                            'group2': g2_name,
                            'p_value': p_val,
                            'p_adj': p_adj,
                            'significant': p_adj < p_threshold
                        })
                        
    except Exception as e:
        print(f"경고: {m_z_bin} 통계 분석 중 오류 - {e}")
        traceback.print_exc() 
        pass

    return main_test_result, posthoc_results


def plot_single_bin(ax, data_bin, m_z_bin, stats_df_bin):
    """
    (수정됨)
    하나의 m/z bin에 대한 개별 플롯(bar + strip)을 생성합니다.
    Scatter plot의 점들이 막대 중앙에 오도록 수정되었습니다.
    """
    
    # 1. Bar plot (막대 그래프)
    sns.barplot(
        data=data_bin,
        x='group',
        y='intensity',
        hue='group',
        palette='Paired',
        errorbar=None,
        ax=ax,
        legend=False 
    )
    
    # 2. Strip plot (산점도) - (수정됨)
    sns.stripplot(
        data=data_bin,
        x='group',
        y='intensity',
        # (수정) hue와 dodge를 제거하여 barplot의 x위치와 일치시킴
        color='black',   # 모든 점을 검은색으로
        legend=False,
        s=5, 
        ax=ax,
        jitter=False # 점들을 중앙에 정렬
    )
    
    # 3. 플롯 제목 및 레이블 설정
    ax.set_title(f"m/z bin: {m_z_bin}", fontsize=10)
    ax.set_xlabel('')
    ax.set_ylabel('Intensity')
    ax.tick_params(axis='x', rotation=45)
    
    # 4. 통계 어노테이션 추가
    significant_pairs = stats_df_bin[stats_df_bin['significant'] == True]
    
    if not significant_pairs.empty:
        max_intensity = data_bin['intensity'].max()
        y_offset = max_intensity * 0.05 
        y = max_intensity + y_offset
        
        x_pos = (len(ax.get_xticklabels()) - 1) / 2.0
        ax.text(x_pos, y, "*", ha='center', va='bottom', fontsize=14, color='red')

def export_to_prism(df_long_roi, roi, config, m_z_bins, output_dir):
    """
    (수정됨)
    GraphPad Prism의 'Grouped' 형식에 맞는 CSV 파일로 데이터를 변환하여 저장합니다.
    유연한 config 구조를 지원합니다.
    """
    try:
        print(f"  -> Prism용 CSV 파일 생성 중: {roi}")
        
        # (수정) 1. 설정에서 유연한 그룹 정보 가져오기
        groups_info_list = config['group_info'] # 딕셔너리의 리스트
        
        # (수정) 3. 컬럼 헤더 생성 (유연하게)
        prism_columns = []
        for group_dict in groups_info_list:
            g = group_dict['name']
            n_range = range(1, group_dict['n_per_group'] + 1)
            for n in n_range:
                prism_columns.append(f"{g}-{n}")
        
        # 2. 최종 Prism DataFrame을 위한 데이터 준비
        prism_data = {} 
        
        # 4. m/z bin을 기준으로 데이터 재정렬
        for m_z_bin in m_z_bins:
            row_data = []
            
            # (수정) 5. 컬럼 생성 순서와 동일하게 데이터 추출
            for group_dict in groups_info_list:
                g = group_dict['name']
                n_range = range(1, group_dict['n_per_group'] + 1)
                for n in n_range:
                    val_series = df_long_roi[
                        (df_long_roi['m_z_bin'] == m_z_bin) &
                        (df_long_roi['group'] == g) &
                        (df_long_roi['n'] == n)
                    ]['intensity']
                    
                    if not val_series.empty:
                        row_data.append(val_series.iloc[0])
                    else:
                        row_data.append(None) 
            
            prism_data[m_z_bin] = row_data

        # 5. DataFrame 생성 및 m/z bin을 인덱스로 설정
        df_prism = pd.DataFrame.from_dict(prism_data, orient='index', columns=prism_columns)
        
        # 6. CSV 파일로 저장
        prism_file_path = os.path.join(output_dir, f"aggregated_data_roi_{roi}_prism.csv")
        df_prism.to_csv(prism_file_path, float_format='%.4f', index=True, index_label='m_z_bin')
        
        print(f"  -> Prism 파일 저장 완료: {prism_file_path}")
    
    except KeyError as e:
        print(f"오류: Prism 변환 중 config.yaml 키 오류 - {e}")
    except Exception as e:
        print(f"오류: Prism용 CSV 파일 생성 중 오류 발생 - {e}")
        traceback.print_exc()

def generate_single_plot(data_bin, roi, m_z_bin, stats_df_bin, group_color_map, groups, output_dir):
    """
    (추가됨) 하나의 m/z bin에 대한 개별 플롯을 생성하고 저장합니다.
    """
    try:
        fig_single, ax_single = plt.subplots(figsize=(6, 5))
        
        plot_single_bin(ax_single, data_bin, m_z_bin, stats_df_bin)
        
        # 개별 플롯에 범례 추가
        handles = [plt.Rectangle((0,0),1,1, color=group_color_map[group]) for group in groups]
        ax_single.legend(handles, groups, title="Groups", bbox_to_anchor=(1.05, 1), loc='upper left')

        bin_filename_safe = str(m_z_bin).replace('.', '_')
        single_plot_path = os.path.join(output_dir, f"plot_roi_{roi}_bin_{bin_filename_safe}.png")
        
        fig_single.savefig(single_plot_path, dpi=150, bbox_inches='tight')
        plt.close(fig_single) 
        
    except Exception as e_single:
        print(f"경고: {m_z_bin} 개별 플롯 저장 중 오류: {e_single}")
        plt.close(fig_single)

def generate_montage_plot(df_long_roi, roi, value_vars, df_posthoc_roi, group_color_map, groups, output_dir):
    """
    (추가됨) 지정된 ROI에 대한 몽타주 플롯을 생성하고 저장합니다.
    """
    print(f"  -> 몽타주 플롯 생성 중: {roi}")
    
    num_bins = len(value_vars)
    cols = int(np.ceil(np.sqrt(num_bins)))
    rows = int(np.ceil(num_bins / cols))
    
    # (개선) GridSpec을 사용하여 범례 공간을 명시적으로 분리
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

    # (개선) 몽타주 범례 추가
    legend_ax = fig.add_subplot(gs[:, -1])
    legend_ax.axis('off')
    handles = [plt.Rectangle((0,0),1,1, color=group_color_map[group]) for group in groups]
    legend_ax.legend(handles, groups, title="Groups", loc='center')

    fig.suptitle(f'MSI Intensity Analysis (ROI: {roi})', fontsize=16, y=1.02)
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    
    montage_plot_path = os.path.join(output_dir, f"plot_montage_roi_{roi}.png")
    fig.savefig(montage_plot_path, dpi=150, bbox_inches='tight')
    print(f"  -> 몽타주 플롯 저장 완료: {montage_plot_path}")
    plt.close(fig)

def generate_html_report(config, output_dir, df_main_stats, df_posthoc_stats):
    """
    분석 결과를 요약하는 HTML 보고서를 생성하고 저장합니다.
    """
    print("\n--- HTML 보고서 생성 시작 ---")
    report_title = config.get('report_settings', {}).get('title', 'Mass Imaging Analysis Report')
    p_threshold = config['statistics_settings']['p_value_threshold']
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{report_title}</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                line-height: 1.6;
                background-color: #f4f7f6;
                color: #333;
                margin: 0;
                padding: 20px;
            }}
            .container {{
                max-width: 1200px;
                margin: 20px auto;
                padding: 30px;
                background-color: #ffffff;
                box-shadow: 0 4px 15px rgba(0,0,0,0.08);
                border-radius: 10px;
            }}
            h1, h2, h3 {{
                color: #2c3e50;
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
                margin-top: 40px;
            }}
            h1 {{ text-align: center; border-bottom: none; }}
            h2.collapsible {{ cursor: pointer; }}
            h2.collapsible::after {{ content: ' ▼'; font-size: 0.8em; color: #3498db; }}
            h2.collapsible.active::after {{ content: ' ▲'; }}
            .content {{ padding: 0 18px; max-height: 0; overflow: hidden; transition: max-height 0.3s ease-out; background-color: #f9f9f9; border-radius: 5px; margin-top: 10px; }}
            .card {{
                background: #fff;
                border: 1px solid #e1e8ed;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 25px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.05);
            }}
            table.dataframe {{
                width: auto;
                max-width: 100%;
                border-collapse: collapse;
                margin: 20px auto;
                font-size: 0.9em;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }}
            table.dataframe th, table.dataframe td {{
                border: 1px solid #ddd;
                padding: 12px 15px;
                text-align: left;
            }}
            table.dataframe thead th {{
                background-color: #3498db;
                color: #ffffff;
                font-weight: bold;
            }}
            table.dataframe tbody tr:nth-of-type(even) {{ background-color: #f3f3f3; }}
            .significant {{ background-color: #ffe5e5 !important; }}
            .no-results {{ color: #777; font-style: italic; padding: 15px; border: 1px dashed #ccc; background-color: #fafafa; text-align: center; }}
            .plot-container {{ text-align: center; margin-top: 20px; }}
            .plot-container img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>{report_title}</h1>
            <p style="text-align: center; color: #777;"><strong>Report generated on:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

            <h2>1. Experiment Configuration Summary</h2>
            <div class="card">
    """

    # 1.1 Group Information
    html_content += "<h3>1.1 Group Information</h3>"
    group_info_df = pd.DataFrame(config['group_info'])
    html_content += group_info_df.to_html(index=False, classes='dataframe', border=0)

    # 1.2 ROI Information
    html_content += "<h3>1.2 ROI Information</h3>"
    roi_info_df = pd.DataFrame(config['roi_info'])
    html_content += roi_info_df.to_html(index=False, classes='dataframe', border=0)

    # 1.3 Binning Settings
    html_content += "<h3>1.3 Binning Settings</h3>"
    binning_settings_df = pd.DataFrame([config['binning_settings']])
    html_content += binning_settings_df.to_html(index=False, classes='dataframe', border=0)

    # 1.4 Statistics Settings
    html_content += "<h3>1.4 Statistics Settings</h3>"
    stats_settings_df = pd.DataFrame([config['statistics_settings']])
    html_content += stats_settings_df.to_html(index=False, classes='dataframe', border=0)
    html_content += "</div>"

    html_content += "<h2>2. Statistical Analysis Results</h2>"
    html_content += "<div class='card'>"

    # 2.1 Significant Results (p-adjusted < p_threshold)
    html_content += f"<h3>2.1 Significant Post-hoc Results (p-adjusted < {p_threshold})</h3>"
    if not df_posthoc_stats.empty:
        significant_results = df_posthoc_stats[df_posthoc_stats['significant'] == True].copy()
        if not significant_results.empty:
            # 유의미한 결과 행 강조
            def highlight_significant_rows(row):
                return ['background-color: #fff0f0' if row['significant'] else '' for _ in row]
            
            html_content += significant_results.style.apply(highlight_significant_rows, axis=1).to_html(index=False, classes='dataframe')
        else:
            html_content += "<p class='no-results'>No statistically significant post-hoc results found at the specified p-value threshold.</p>"
    else:
        html_content += "<p class='no-results'>No post-hoc statistical results were generated.</p>"

    # 2.2 Full Main Statistical Results
    html_content += "<h3>2.2 Full Main Statistical Results</h3>"
    html_content += df_main_stats.to_html(index=False, classes='dataframe')

    # 2.3 Full Post-hoc Statistical Results
    html_content += "<h3>2.3 Full Post-hoc Statistical Results</h3>"
    if not df_posthoc_stats.empty:
        html_content += df_posthoc_stats.to_html(index=False, classes='dataframe')
    else:
        html_content += "<p class='no-results'>No post-hoc statistical results were generated.</p>"

    html_content += "<h2>3. Visualizations</h2>"

    # 3.1 Montage Plots
    if 'roi_info' in config:
        for roi_dict in config['roi_info']:
            roi_name = roi_dict['name']
            plot_filename = f"plot_montage_roi_{roi_name}.png"
            # HTML 보고서와 이미지가 같은 디렉토리에 저장되므로 파일명만으로 참조 가능
            if os.path.exists(os.path.join(output_dir, plot_filename)):
                html_content += f"""
                <h3>3.1 Montage Plot for ROI: {roi_name}</h3>
                <div class="plot-container">
                    <img src="{plot_filename}" alt="Montage Plot for {roi_name}">
                </div>
                """
            else:
                html_content += f"<p class='no-results'>Montage plot for ROI '{roi_name}' not found.</p>"
    else:
        html_content += "<p class='no-results'>No ROI information found in config to generate plots.</p>"

    html_content += """
    </body>
    </html>
    """

    report_filepath = os.path.join(output_dir, "analysis_report.html")
    with open(report_filepath, 'w', encoding='utf-8') as f:
        f.write(html_content)
    print(f"HTML report generated: {report_filepath}")


def main(config_path='config/config.yaml'):
    """
    (수정됨)
    통계 분석 및 시각화 메인 함수
    유연한 config 구조를 지원합니다.
    config_path를 인자로 받아 해당 설정으로 분석을 수행합니다.
    """
    print("--- 1. 설정 파일 로드 ---")
    config = parsing.load_yaml(config_path)
    if config is None:
        return 

    try:
        output_dir = config['output_dir']
        # (수정) 딕셔너리 리스트에서 'name'만 추출
        rois = [r['name'] for r in config['roi_info']]
        stats_settings = config['statistics_settings']
        test_type = stats_settings['test_type']
        p_threshold = stats_settings['p_value_threshold']
        # (수정) 딕셔너리 리스트에서 'name'만 추출 (이것이 그룹 순서를 결정)
        groups = [g['name'] for g in config['group_info']]
        
    except KeyError as e:
        print(f"오류: config.yaml 파일에 필요한 키가 없습니다: {e}")
        return
        
    print("--- 2. 집계 데이터 로드 ---")
    agg_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    if not os.path.exists(agg_csv_path):
        print(f"오류: 집계 파일 '{agg_csv_path}'을(를) 찾을 수 없습니다.")
        print("'aggregate.py'가 먼저 실행되었는지 확인하세요.")
        return
        
    df_agg = pd.read_csv(agg_csv_path)
    
    # 3. 데이터 "Melt" (Wide -> Long format)
    id_vars = ['group', 'n', 'roi']
    value_vars = [col for col in df_agg.columns if col not in id_vars]
    
    df_long = df_agg.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        var_name='m_z_bin',
        value_name='intensity'
    )
    
    try:
        # (수정) config.yaml에서 읽어온 'groups' 리스트를 순서로 지정
        df_long['group'] = pd.Categorical(df_long['group'], categories=groups, ordered=True)
    except Exception as e:
        print(f"경고: 그룹 순서 적용 중 오류: {e}")
    
    print(f"데이터를 Long format으로 변환 (총 {len(df_long)} 행)")
    
    # 4. ROI 및 m/z bin별 통계 분석 및 플로팅
    all_main_stats = []
    all_posthoc_stats = []
    
    # (수정) config.yaml 순서에 맞게 색상 매핑
    palette = sns.color_palette("Paired", n_colors=len(groups))
    group_color_map = dict(zip(groups, palette))
    
    print(f"\n--- 3. ROI {len(rois)}개에 대한 분석 시작 ---")
    for roi in rois:
        print(f"\n분석 중인 ROI: {roi}")
        
        df_long_roi = df_long[df_long['roi'] == roi].copy()
        if df_long_roi.empty:
            print(f"  -> 경고: {roi}에 대한 데이터가 없습니다.")
            continue
        
        # 4-1. 통계 분석 (m/z bin별로 수행)
        for m_z_bin in value_vars: 
            main_result, posthoc_results = run_statistics(
                df_long_roi, m_z_bin, roi, test_type, p_threshold
            )
            all_main_stats.append(main_result)
            all_posthoc_stats.extend(posthoc_results) 

        df_posthoc_roi = pd.DataFrame(all_posthoc_stats) 
        
        # --- 몽타주 플롯 생성 로직 ---
        print(f"  -> 몽타주 플롯 생성 중: {roi}")
        
        num_bins = len(value_vars)
        cols = int(np.ceil(np.sqrt(num_bins)))
        rows = int(np.ceil(num_bins / cols))
        
        fig = plt.figure(figsize=(cols * 5, rows * 4))
        gs = gridspec.GridSpec(rows, cols + 1, figure=fig, width_ratios=[1] * cols + [0.5])
        
        plot_axes = [] 
        for r in range(rows):
            for c in range(cols):
                idx = r * cols + c
                if idx < num_bins:
                    plot_axes.append(fig.add_subplot(gs[r, c]))
        
        for idx, (m_z_bin, ax) in enumerate(zip(value_vars, plot_axes)):
            data_bin = df_long_roi[df_long_roi['m_z_bin'] == m_z_bin]
            
            stats_df_bin = pd.DataFrame() 
            if not df_posthoc_roi.empty: 
                stats_df_bin = df_posthoc_roi[
                    (df_posthoc_roi['m_z_bin'] == m_z_bin) & 
                    (df_posthoc_roi['roi'] == roi)
                ] if not df_posthoc_roi.empty else pd.DataFrame()

            plot_single_bin(ax, data_bin, m_z_bin, stats_df_bin)
            generate_single_plot(data_bin, roi, m_z_bin, stats_df_bin, group_color_map, groups, output_dir)

        # (개선) 몽타주 플롯 생성 로직을 별도 함수로 분리
        generate_montage_plot(df_long_roi, roi, value_vars, df_posthoc_roi, group_color_map, groups, output_dir)
        
        # 4-4. Prism용 CSV 내보내기
        export_to_prism(df_long_roi, roi, config, value_vars, output_dir)

    # 5. 통계 결과 CSV로 저장
    print("\n--- 4. 통계 결과 저장 ---")
    try:
        df_main_stats = pd.DataFrame(all_main_stats)
        main_stats_path = os.path.join(output_dir, "statistical_results_main.csv")
        df_main_stats.to_csv(main_stats_path, index=False, float_format='%.4e')
        print(f"메인 통계 결과 저장 완료: {main_stats_path}")

        if all_posthoc_stats: 
            df_posthoc_stats = pd.DataFrame(all_posthoc_stats)
            posthoc_stats_path = os.path.join(output_dir, "statistical_results_posthoc.csv")
            df_posthoc_stats.to_csv(posthoc_stats_path, index=False, float_format='%.4e')
            print(f"Post-hoc 통계 결과 저장 완료: {posthoc_stats_path}")
        else:
            # (수정) Post-hoc 결과가 없을 때도 빈 DataFrame을 생성하여 NameError 방지
            df_posthoc_stats = pd.DataFrame()
            print("Post-hoc 통계 결과가 없습니다 (유의미한 차이 없음).")

    except Exception as e:
        print(f"오류: 통계 결과 파일 저장 중 오류 발생 - {e}")

    # 6. HTML 보고서 생성
    generate_html_report(config, output_dir, df_main_stats, df_posthoc_stats)

    print("\n--- 모든 분석 및 시각화가 완료되었습니다. ---")

if __name__ == '__main__':
    # 스크립트를 직접 실행할 경우 기본 config 파일을 사용
    main(config_path='config/config.yaml')