import os
import pandas as pd
import yaml
import itertools
import traceback

# 통계 및 시각화 라이브러리
import seaborn as sns
import matplotlib.pyplot as plt
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
        # sys.exit(1) 대신 return을 사용 (Jupyter 호환)
        # main() 함수가 아니므로 이 파일 임포트 시 오류 발생 가능
        pass 

CONFIG_FILE = 'config/config.yaml'

def run_statistics(data_long_roi, m_z_bin, roi, test_type, p_threshold):
    """
    특정 ROI와 m/z bin에 대해 통계 분석을 수행합니다.
    (실제 통계 로직으로 완전 대체)
    """
    
    # 1. 분석할 데이터 필터링
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
    posthoc_results = [] # 사후 검정 결과를 담을 리스트

    if n_groups < 2:
        print(f"  -> 경고: {m_z_bin}에 그룹이 1개뿐이라 통계 불가.")
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
                stat, p_val = ttest_ind(g1_data, g2_data, equal_var=False) # Welch's T-test
            else:
                test_name = 'Mann-Whitney U'
                stat, p_val = mannwhitneyu(g1_data, g2_data, alternative='two-sided')
            
            main_test_result['test_name'] = test_name
            main_test_result['p_value'] = p_val
            main_test_result['significant'] = p_val < p_threshold
            
            # 2개 그룹 비교는 main test 자체가 post-hoc임
            posthoc_results.append({
                'roi': roi,
                'm_z_bin': m_z_bin,
                'test_name': test_name,
                'group1': g1_name,
                'group2': g2_name,
                'p_value': p_val,
                'p_adj': p_val, # p-value가 1개이므로 조정 없음
                'significant': p_val < p_threshold
            })

        else:
            # --- 3개 이상 그룹 비교 (ANOVA or Kruskal-Wallis) ---
            if test_type == 'parametric':
                # ANOVA
                test_name = 'ANOVA'
                stat, p_val = f_oneway(*all_group_data)
            else:
                # Kruskal-Wallis
                test_name = 'Kruskal-Wallis'
                stat, p_val = kruskal(*all_group_data)
            
            main_test_result['test_name'] = test_name
            main_test_result['p_value'] = p_val
            main_test_result['significant'] = p_val < p_threshold

            # --- Post-hoc (사후 검정) ---
            if main_test_result['significant']:
                if test_type == 'parametric':
                    # Tukey's HSD
                    # Tukey는 long-form 데이터가 필요
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
                            'p_value': row['p-adj'], # Tukey는 p-adj를 p-value로 사용
                            'p_adj': row['p-adj'],
                            'significant': row['reject'] # True/False
                        })
                else:
                    # Non-parametric post-hoc: Pairwise Mann-Whitney U with Bonferroni correction
                    posthoc_test_name = 'Mann-Whitney (Bonferroni)'
                    all_pairs = list(itertools.combinations(groups, 2))
                    num_comparisons = len(all_pairs)
                    
                    for g1_name, g2_name in all_pairs:
                        g1_data = data[data['group'] == g1_name]['intensity']
                        g2_data = data[data['group'] == g2_name]['intensity']
                        
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

    return main_test_result, posthoc_results


def create_plots(data_long_roi, roi, stats_df, p_threshold, output_dir):
    """
    통계 결과(post-hoc)를 포함하여 ROI별로 플롯을 생성하고 저장합니다.
    (Scatter plot이 그룹별로 분리(dodge)되도록 수정)
    """
    try:
        print(f"  -> 플롯 생성 중: {roi}")
        
        # 1. 기본 플롯 (catplot): 막대와 기본 범례 생성
        g = sns.catplot(
            data=data_long_roi,
            x='m_z_bin',
            y='intensity',
            hue='group',
            col='roi',
            kind='bar',
            errorbar=None,
            palette='Paired',
            legend=True,
            aspect=2
        )
        g.despine(left=True)
        
        # 2. 오버레이 플롯 (stripplot): 개별 데이터 포인트를 검은색 점으로 추가
        g.map_dataframe(
            sns.stripplot,
            x='m_z_bin',
            y='intensity',
            # (핵심 수정)
            hue='group',     # 'dodge'가 그룹을 인식하도록 hue 전달
            dodge=True,      # 그룹별로 점 분리
            palette='dark:black',   # 점 색상을 검은색으로 통일
            legend=False,    # 이 플롯은 범례에 항목 추가 안 함
            s=5,             # 점 크기
        )

        # 3. 플롯 제목 및 축 레이블 설정
        g.fig.suptitle(f'Intensity by Group and m/z bin (ROI: {roi})', y=1.03)
        g.set_axis_labels("m/z Bin", "Mean Intensity")
        g.set_xticklabels(rotation=45)

        # 4. 통계적 유의성 표시 (p-value 어노테이션)
        try:
            # 유의미한 쌍(pair)만 필터링
            # (stats_df는 df_posthoc임)
            significant_pairs = stats_df[
                (stats_df['significant'] == True) & 
                (stats_df['roi'] == roi)
            ]
            
            if not significant_pairs.empty:
                print(f"  -> {roi}: {len(significant_pairs)}개의 유의미한 쌍에 대한 어노테이션 추가 시도")
                
                m_z_bins_order = data_long_roi['m_z_bin'].unique()
                ax = g.axes[0, 0] # (col='roi'가 하나라고 가정)
                
                for i, m_z_bin in enumerate(m_z_bins_order):
                    pairs_in_bin = significant_pairs[significant_pairs['m_z_bin'] == m_z_bin]
                    if pairs_in_bin.empty:
                        continue

                    # 해당 m/z bin의 intensity 최대값 (y축 위치)
                    max_intensity = data_long_roi[data_long_roi['m_z_bin'] == m_z_bin]['intensity'].max()
                    y_offset = max_intensity * 0.05 # 5% 높이
                    y = max_intensity + y_offset

                    # (간단하게 '*'로 유의미함만 표시)
                    if not pairs_in_bin.empty:
                         ax.text(i, y, "*", ha='center', va='bottom', fontsize=14, color='red')

        except KeyError as e:
            # 이 오류가 발생했던 지점
            print(f"경고: 플롯 어노테이션 중 KeyError - {e}. stats_df 컬럼 확인 필요.")
            print(f"  stats_df columns: {stats_df.columns.tolist() if not stats_df.empty else 'DataFrame is empty'}")
            
        except Exception as e_stats:
            print(f"경고: 플롯에 통계 어노테이션 추가 중 오류 발생 - {e_stats}")
            traceback.print_exc()

        # 5. 플롯 저장
        plot_path = os.path.join(output_dir, f"plot_roi_{roi}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"  -> 플롯 저장 완료: {plot_path}")
        plt.close(g.fig) # 메모리 해제

    except Exception as e:
        print(f"오류: 플롯 생성 중 심각한 오류 발생 - {e}")
        traceback.print_exc() # 디버깅을 위한 상세 오류 출력


def export_to_prism(df_long_roi, roi, config, m_z_bins, output_dir):
    """
    GraphPad Prism의 'Grouped' 형식에 맞는 CSV 파일로 데이터를 변환하여 저장합니다.
    (이전 버전과 동일, 변경 없음)
    """
    try:
        print(f"  -> Prism용 CSV 파일 생성 중: {roi}")
        
        # 1. 설정에서 그룹 정보 가져오기
        groups = config['group_info']['name_of_group']
        n_count = config['n_per_group']
        n_range = range(1, n_count + 1)
        
        # 2. 최종 Prism DataFrame을 위한 데이터 준비
        prism_data = {} # Key: m/z bin, Value: [list of intensities]
        
        # 3. 컬럼 헤더 생성 (예: 'acetate-1', 'acetate-2', ...)
        prism_columns = []
        for g in groups:
            for n in n_range:
                prism_columns.append(f"{g}-{n}")
        
        # 4. m/z bin을 기준으로 데이터 재정렬
        for m_z_bin in m_z_bins:
            row_data = []
            
            # config에 정의된 그룹 순서대로 데이터 추출
            for g in groups:
                for n in n_range:
                    # df_long_roi에서 일치하는 데이터 찾기
                    val_series = df_long_roi[
                        (df_long_roi['m_z_bin'] == m_z_bin) &
                        (df_long_roi['group'] == g) &
                        (df_long_roi['n'] == n)
                    ]['intensity']
                    
                    if not val_series.empty:
                        row_data.append(val_series.iloc[0])
                    else:
                        row_data.append(None) # 데이터가 없는 경우
            
            prism_data[m_z_bin] = row_data

        # 5. DataFrame 생성 및 m/z bin을 인덱스로 설정
        df_prism = pd.DataFrame.from_dict(prism_data, orient='index', columns=prism_columns)
        
        # 6. CSV 파일로 저장
        # Prism이 인식하도록 인덱스(m/z bin)도 저장 (index=True, index_label='m_z_bin')
        prism_file_path = os.path.join(output_dir, f"aggregated_data_roi_{roi}_prism.csv")
        df_prism.to_csv(prism_file_path, float_format='%.4f', index=True, index_label='m_z_bin')
        
        print(f"  -> Prism 파일 저장 완료: {prism_file_path}")
    
    except KeyError as e:
        print(f"오류: Prism 변환 중 config.yaml 키 오류 - {e}")
    except Exception as e:
        print(f"오류: Prism용 CSV 파일 생성 중 오류 발생 - {e}")
        traceback.print_exc() # 디버깅을 위한 상세 오류 출력


def main():
    """
    통계 분석 및 시각화 메인 함수
    """
    print("--- 1. 설정 파일 로드 ---")
    config = parsing.load_yaml(CONFIG_FILE)
    if config is None:
        return # (Jupyter 호환)

    try:
        output_dir = config['output_dir']
        rois = config['roi_info']['name_of_roi']
        stats_settings = config['statistics_settings']
        test_type = stats_settings['test_type']
        p_threshold = stats_settings['p_value_threshold']
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
    # m/z bin 컬럼(value_vars)을 자동으로 식별
    id_vars = ['group', 'n', 'roi']
    value_vars = [col for col in df_agg.columns if col not in id_vars]
    
    df_long = df_agg.melt(
        id_vars=id_vars,
        value_vars=value_vars,
        var_name='m_z_bin',      # 새로운 m/z bin 컬럼 이름
        value_name='intensity'   # 새로운 intensity 컬럼 이름
    )
    
    print(f"데이터를 Long format으로 변환 (총 {len(df_long)} 행)")
    
    # 4. ROI 및 m/z bin별 통계 분석 및 플로팅
    all_main_stats = []
    all_posthoc_stats = []
    
    print(f"\n--- 3. ROI {len(rois)}개에 대한 분석 시작 ---")
    for roi in rois:
        print(f"\n분석 중인 ROI: {roi}")
        
        # 현재 ROI에 해당하는 데이터만 필터링
        df_long_roi = df_long[df_long['roi'] == roi].copy()
        if df_long_roi.empty:
            print(f"  -> 경고: {roi}에 대한 데이터가 없습니다.")
            continue
        
        # 4-1. 통계 분석 (m/z bin별로 수행)
        for m_z_bin in value_vars: # value_vars가 m/z bin 이름 리스트임
            main_result, posthoc_results = run_statistics(
                df_long_roi, m_z_bin, roi, test_type, p_threshold
            )
            all_main_stats.append(main_result)
            all_posthoc_stats.extend(posthoc_results) # extend로 리스트 병합

        # 4-2. 통계 결과를 DataFrame으로 변환 (플로팅에 사용)
        df_posthoc_roi = pd.DataFrame(all_posthoc_stats)
        
        # 4-3. 플롯 생성 (통계 결과 전달)
        create_plots(df_long_roi, roi, df_posthoc_roi, p_threshold, output_dir)
        
        # 4-4. Prism용 CSV 내보내기
        export_to_prism(df_long_roi, roi, config, value_vars, output_dir)

    # 5. 통계 결과 CSV로 저장
    print("\n--- 4. 통계 결과 저장 ---")
    try:
        df_main_stats = pd.DataFrame(all_main_stats)
        main_stats_path = os.path.join(output_dir, "statistical_results_main.csv")
        df_main_stats.to_csv(main_stats_path, index=False, float_format='%.4e')
        print(f"메인 통계 결과 저장 완료: {main_stats_path}")

        if all_posthoc_stats: # Post-hoc 결과가 있을 때만 저장
            df_posthoc_stats = pd.DataFrame(all_posthoc_stats)
            posthoc_stats_path = os.path.join(output_dir, "statistical_results_posthoc.csv")
            df_posthoc_stats.to_csv(posthoc_stats_path, index=False, float_format='%.4e')
            print(f"Post-hoc 통계 결과 저장 완료: {posthoc_stats_path}")
        else:
            print("Post-hoc 통계 결과가 없습니다 (유의미한 차이 없음).")

    except Exception as e:
        print(f"오류: 통계 결과 파일 저장 중 오류 발생 - {e}")

    print("\n--- 모든 분석 및 시각화가 완료되었습니다. ---")

if __name__ == '__main__':
    # 이 스크립트를 직접 실행할 경우
    # (Jupyter Notebook에서 import하는 것이 기본)
    main()