import os
import pandas as pd
import yaml
import itertools
import numpy as np

# 'parsing.py'의 load_yaml 함수를 재사용합니다.
try:
    import parsing
except ImportError:
    print("오류: 'parsing.py' 파일을 찾을 수 없습니다.")
    print("aggregate.py와 parsing.py가 동일한 디렉토리에 있는지 확인하세요.")
    exit(1)

CONFIG_FILE = 'config/config.yaml'

def main():
    """
    '_mean_intensities.csv' 파일들을 집계(aggregate)하여
    num_serial (technical replicates) 간의 평균을 계산하고
    하나의 최종 CSV 파일로 저장합니다.
    """
    print("--- 1. 설정 파일 로드 ---")
    config = parsing.load_yaml(CONFIG_FILE)
    if config is None:
        return

    try:
        # 1-1. 설정 정보 가져오기
        output_dir = config['output_dir']
        groups = config['group_info']['name_of_group']
        n_range = range(1, config['n_per_group'] + 1)
        s_range = range(1, config['num_serial'] + 1)
        rois = config['roi_info']['name_of_roi']
        
    except KeyError as e:
        print(f"오류: config.yaml 파일에 필요한 키가 없습니다: {e}")
        return
        
    print(f"'{output_dir}' 디렉토리에서 '_mean_intensities.csv' 파일들을 집계합니다.")

    # 1-2. 집계의 기준이 될 조합 생성 (num_serial 제외)
    # (group, n, roi)의 12개 조합 (4 * 3 * 1)
    aggregation_groups = list(itertools.product(groups, n_range, rois))
    
    all_aggregated_data = []
    m_z_columns = None # m/z 빈 컬럼 이름을 저장할 변수

    print(f"\n--- 2. 총 {len(aggregation_groups)}개 그룹 집계 시작 ---")

    # 2. (group, n, roi) 조합을 기준으로 루프 (총 12번)
    for (group, n, roi) in aggregation_groups:
        
        serial_data_to_average = [] # s=1, s=2 등의 데이터를 담을 리스트
        
        # 3. 'num_serial' (s_range)을 기준으로 파일 검색
        for s in s_range:
            # main.py에서 생성한 파일명 템플릿과 동일해야 함
            base_imzml_name = f"{group} {n}-{s} {roi}-total ion count"
            mean_csv_path = os.path.join(output_dir, f"{base_imzml_name}_mean_intensities.csv")
            
            if os.path.exists(mean_csv_path):
                try:
                    df_mean = pd.read_csv(mean_csv_path)
                    serial_data_to_average.append(df_mean)
                    
                    # m/z 컬럼 이름은 첫 번째 파일에서 한 번만 가져옵니다.
                    if m_z_columns is None:
                        m_z_columns = df_mean.columns.tolist()
                        print(f"m/z 빈(컬럼) {len(m_z_columns)}개 확인: {m_z_columns}")
                        
                except Exception as e:
                    print(f"오류: {mean_csv_path} 파일 읽기 실패 - {e}")
            else:
                print(f"경고: 집계에 필요한 파일을 찾을 수 없음: {mean_csv_path}")

        # 4. 'num_serial' 간 평균 계산
        if not serial_data_to_average:
            print(f"  -> 처리 실패: {group} {n} {roi} 그룹에 대한 데이터 없음.")
            continue
            
        # s=1, s=2에 해당하는 DataFrame들을 하나로 합침 (세로로)
        df_concat = pd.concat(serial_data_to_average)
        
        # 합쳐진 DataFrame의 m/z 컬럼별 평균을 계산 (axis=0)
        df_averaged = df_concat.mean(axis=0) # Pandas Series 반환
        
        # 5. 결과 저장 (Series를 dictionary로 변환)
        agg_row_data = df_averaged.to_dict()
        
        # 식별자(Identifier) 컬럼 추가
        agg_row_data['group'] = group
        agg_row_data['n'] = n
        agg_row_data['roi'] = roi
        
        all_aggregated_data.append(agg_row_data)
        print(f"  -> 처리 완료: {group} {n} {roi} (파일 {len(serial_data_to_average)}개 평균)")

    # 6. 최종 DataFrame 생성 및 저장
    if not all_aggregated_data:
        print("오류: 집계할 데이터가 없습니다. 'main.py'가 정상적으로 실행되었는지 확인하세요.")
        return
        
    df_final = pd.DataFrame(all_aggregated_data)
    
    # 컬럼 순서 재정렬 (식별자 컬럼을 맨 앞으로)
    if m_z_columns:
        identifier_cols = ['group', 'n', 'roi']
        final_cols = identifier_cols + m_z_columns
        # 혹시 모를 불일치 방지
        final_cols = [col for col in final_cols if col in df_final.columns] 
        df_final = df_final[final_cols]
    else:
        print("경고: m/z 컬럼이 설정되지 않아 컬럼 순서를 재정렬할 수 없습니다.")

    # 7. 최종 CSV 파일로 저장
    output_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    df_final.to_csv(output_csv_path, index=False, float_format='%.4f')
    
    print(f"\n--- 집계가 완료되었습니다. ---")
    print(f"최종 집계 데이터 (형태: {df_final.shape})가 다음 파일에 저장되었습니다:")
    print(f"{output_csv_path}")

if __name__ == '__main__':
    main()
