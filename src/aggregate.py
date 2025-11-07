import os
import pandas as pd
import yaml
import itertools
import numpy as np

# 'src' 폴더의 'parsing.py'의 load_yaml 함수를 재사용합니다.
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

def main(config_path='config/config.yaml'):
    """
    (수정됨)
    '_mean_intensities.csv' 파일들을 집계(aggregate)하여
    유연한 config 설정에 따라 num_serial 간의 평균을 계산하고
    하나의 최종 CSV 파일로 저장합니다.
    config_path를 인자로 받아 해당 설정으로 집계를 수행합니다.
    """
    print("--- 1. 설정 파일 로드 ---")
    config = parsing.load_yaml(config_path)
    if config is None:
        return # (Jupyter 호환)

    try:
        # 1-1. (수정) 유연한 설정 정보 가져오기
        output_dir = config['output_dir']
        groups_info_list = config['group_info'] # 딕셔너리의 리스트
        rois_info_list = config['roi_info']   # 딕셔너리의 리스트
        
    except KeyError as e:
        print(f"오류: config.yaml 파일에 필요한 키가 없습니다: {e}")
        return # (Jupyter 호환)
        
    print(f"'{output_dir}' 디렉토리에서 '_mean_intensities.csv' 파일들을 집계합니다.")

    all_aggregated_data = []
    m_z_columns = None # m/z 빈 컬럼 이름을 저장할 변수

    # 1-2. (수정) 집계 기준이 될 조합을 config 리스트 순회를 통해 생성
    print(f"\n--- 2. 총 {len(groups_info_list) * sum(g['n_per_group'] for g in groups_info_list) * len(rois_info_list)}개 조합(group*n*roi) 집계 시도 ---")

    # 2. (수정) 'group_info' 리스트를 순회
    for group_dict in groups_info_list:
        group_name = group_dict['name']
        n_range = range(1, group_dict['n_per_group'] + 1)
        s_range = range(1, group_dict['num_serial'] + 1)
        
        # 'roi_info' 리스트를 순회
        for roi_dict in rois_info_list:
            roi_name = roi_dict['name']
            
            # 'n_per_group' (biological replicates)을 순회
            for n in n_range:
                
                serial_data_to_average = [] # s=1, s=2 등의 데이터를 담을 리스트
                
                # 3. 'num_serial' (technical replicates)을 기준으로 파일 검색
                for s in s_range:
                    # main.py에서 생성한 파일명 템플릿과 동일
                    base_imzml_name = f"{group_name} {n}-{s} {roi_name}-total ion count"
                    mean_csv_path = os.path.join(output_dir, f"{base_imzml_name}_mean_intensities.csv")
                    
                    if os.path.exists(mean_csv_path):
                        try:
                            df_mean = pd.read_csv(mean_csv_path)
                            serial_data_to_average.append(df_mean)
                            
                            if m_z_columns is None:
                                m_z_columns = df_mean.columns.tolist()
                                print(f"m/z 빈(컬럼) {len(m_z_columns)}개 확인: {m_z_columns[0]}...")
                                
                        except Exception as e:
                            print(f"오류: {mean_csv_path} 파일 읽기 실패 - {e}")
                    else:
                        print(f"경고: 집계에 필요한 파일을 찾을 수 없음: {mean_csv_path}")

                # 4. 'num_serial' 간 평균 계산
                if not serial_data_to_average:
                    print(f"  -> 처리 실패: {group_name} {n} {roi_name} 그룹에 대한 데이터 없음.")
                    continue
                    
                df_concat = pd.concat(serial_data_to_average)
                df_averaged = df_concat.mean(axis=0) # Pandas Series 반환
                
                # 5. 결과 저장 (Series를 dictionary로 변환)
                agg_row_data = df_averaged.to_dict()
                
                # 식별자(Identifier) 컬럼 추가
                agg_row_data['group'] = group_name
                agg_row_data['n'] = n
                agg_row_data['roi'] = roi_name
                
                all_aggregated_data.append(agg_row_data)
                print(f"  -> 처리 완료: {group_name} {n} {roi_name} (파일 {len(serial_data_to_average)}개 평균)")

    # 6. 최종 DataFrame 생성 및 저장
    if not all_aggregated_data:
        print("오류: 집계할 데이터가 없습니다. 'main.py'가 정상적으로 실행되었는지 확인하세요.")
        return # (Jupyter 호환)
        
    df_final = pd.DataFrame(all_aggregated_data)
    
    # 컬럼 순서 재정렬 (식별자 컬럼을 맨 앞으로)
    if m_z_columns:
        identifier_cols = ['group', 'n', 'roi']
        final_cols = identifier_cols + m_z_columns
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
    # 스크립트를 직접 실행할 경우 기본 config 파일을 사용
    main(config_path='config/config.yaml')