import os
import pandas as pd
import yaml
import itertools
import numpy as np

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

def main(config_path='config/config.yaml'):
    """데이터 집계 실행 함수"""
    logger.info("========== Step 2: 데이터 집계 시작 ==========")
    
    config = parsing.load_yaml(config_path)
    if config is None:
        return

    try:
        output_dir = config['output_dir']
        groups_info_list = config['group_info']
        rois_info_list = config['roi_info']
    except KeyError as e:
        logger.error(f"config.yaml 키 오류: {e}")
        return
        
    logger.info(f"'{output_dir}' 폴더의 데이터를 집계합니다.")

    all_aggregated_data = []
    m_z_columns = None

    # 예상 조합 수 계산
    total_combinations = len(groups_info_list) * sum(g['n_per_group'] for g in groups_info_list) * len(rois_info_list)
    # logger.info(f"총 {total_combinations}개 조합(group*n*roi) 집계 시도")

    for group_dict in groups_info_list:
        group_name = group_dict['name']
        n_range = range(1, group_dict['n_per_group'] + 1)
        s_range = range(1, group_dict['num_serial'] + 1)
        
        for roi_dict in rois_info_list:
            roi_name = roi_dict['name']
            
            for n in n_range:
                serial_data_to_average = []
                
                for s in s_range:
                    base_imzml_name = f"{group_name} {n}-{s} {roi_name}-total ion current"
                    mean_csv_path = os.path.join(output_dir, f"{base_imzml_name}_mean_intensities.csv")
                    
                    if os.path.exists(mean_csv_path):
                        try:
                            df_mean = pd.read_csv(mean_csv_path)
                            serial_data_to_average.append(df_mean)
                            
                            if m_z_columns is None:
                                m_z_columns = df_mean.columns.tolist()
                                # logger.info(f"m/z 컬럼 감지됨: {len(m_z_columns)}개")
                        except Exception as e:
                            logger.error(f"파일 읽기 실패 ({mean_csv_path}): {e}")
                    else:
                        logger.warning(f"파일 누락됨: {mean_csv_path}")

                if not serial_data_to_average:
                    logger.warning(f"데이터 없음 - 건너뜀: {group_name} n={n} {roi_name}")
                    continue
                    
                df_concat = pd.concat(serial_data_to_average)
                df_averaged = df_concat.mean(axis=0)
                
                agg_row_data = df_averaged.to_dict()
                agg_row_data['group'] = group_name
                agg_row_data['n'] = n
                agg_row_data['roi'] = roi_name
                
                all_aggregated_data.append(agg_row_data)
                logger.info(f"집계 완료: {group_name} n={n} {roi_name} ({len(serial_data_to_average)} files)")

    if not all_aggregated_data:
        logger.error("집계된 데이터가 없습니다. Step 1이 정상적으로 수행되었는지 확인하세요.")
        return
        
    df_final = pd.DataFrame(all_aggregated_data)
    
    if m_z_columns:
        identifier_cols = ['group', 'n', 'roi']
        final_cols = identifier_cols + m_z_columns
        final_cols = [col for col in final_cols if col in df_final.columns] 
        df_final = df_final[final_cols]

    output_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    df_final.to_csv(output_csv_path, index=False, float_format='%.4f')
    
    logger.info(f"최종 집계 파일 저장 완료: {output_csv_path}")
    logger.info("========== Step 2 완료 ==========")

if __name__ == '__main__':
    main(config_path='config/config.yaml')