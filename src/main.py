import os
import itertools
import yaml

# (NEW) 로깅 유틸리티 임포트
try:
    import parsing
    from utils.logging_utils import get_logger
except ImportError:
    try:
        from . import parsing
        from .utils.logging_utils import get_logger
    except ImportError:
        print("오류: 모듈 임포트 실패.")
        pass 

CONFIG_FILE = 'config/config.yaml'
logger = get_logger(__name__) # 로거 초기화

def generate_expected_files(config):
    """
    config.yaml의 유연한 '딕셔너리 리스트' 구조를 바탕으로
    예상되는 모든 .imzML 파일의 전체 경로 리스트를 생성합니다.
    """
    try:
        data_dir = config['data_dir']
        groups_info_list = config['group_info']
        rois_info_list = config['roi_info']
        
        expected_files = []
        
        for group_dict in groups_info_list:
            group_name = group_dict['name']
            n_count = group_dict['n_per_group']
            s_count = group_dict['num_serial']
            
            n_range = range(1, n_count + 1)
            s_range = range(1, s_count + 1)
            
            for roi_dict in rois_info_list:
                roi_name = roi_dict['name']
                
                for n, s in itertools.product(n_range, s_range):
                    # 파일명 패턴 (추후 config에서 가져오도록 개선 가능)
                    file_name = f"{group_name} {n}-{s} {roi_name}-total ion current.imzML"
                    full_path = os.path.join(data_dir, file_name)
                    expected_files.append(full_path)
            
        return expected_files

    except KeyError as e:
        logger.error(f"config.yaml 키 오류: {e}")
        return None
    except Exception as e:
        logger.error(f"파일 목록 생성 중 오류: {e}")
        return None

def validate_files(file_list):
    logger.info(f"파일 유효성 검사 시작 (총 {len(file_list)}개 파일 예상)")
    missing_files = []
    for f_path in file_list:
        if not os.path.exists(f_path):
            missing_files.append(f_path)
            
    if missing_files:
        logger.error(f"다음 {len(missing_files)}개의 파일을 찾을 수 없습니다:")
        for missing in missing_files:
            logger.error(f"  - {missing}")
        logger.error(".imzML 파일과 .ibd 파일이 모두 존재하는지 확인하세요.")
        return False
    
    logger.info("모든 예상 파일이 존재합니다.")
    return True

def main(config_path='config/config.yaml'):
    """메인 파이프라인 실행 함수"""
    logger.info("========== Step 1: 데이터 파싱 시작 ==========")
    
    config = parsing.load_yaml(config_path)
    if config is None:
        return

    if 'binning_settings' not in config:
        logger.error("설정 파일에 'binning_settings'가 없습니다.")
        return
        
    master_bins, bin_names = parsing.load_master_bins(config['binning_settings'])
    if master_bins is None:
        return

    file_list = generate_expected_files(config)
    if file_list is None:
        return

    if not validate_files(file_list):
        return
        
    output_dir = config.get('output_dir', 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"총 {len(file_list)}개 파일 처리를 시작합니다.")
    
    for i, imzml_file in enumerate(file_list):
        logger.info(f"--- 파일 처리 [{i+1}/{len(file_list)}] ---")
        parsing.process_imzml_with_bins(imzml_file, master_bins, bin_names, output_dir)
        
    logger.info("========== Step 1 완료: 모든 파일 처리됨 ==========")

if __name__ == '__main__':
    main(config_path='config/config.yaml')