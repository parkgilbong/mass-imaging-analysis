import os
import itertools
import yaml

# 'src' 폴더에 함께 있는 parsing.py에서 필요한 함수들을 import합니다.
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

def generate_expected_files(config):
    """
    (수정됨)
    config.yaml의 유연한 '딕셔너리 리스트' 구조를 바탕으로
    예상되는 모든 .imzML 파일의 전체 경로 리스트를 생성합니다.
    """
    try:
        data_dir = config['data_dir']
        
        # (수정) group_info와 roi_info를 리스트로 직접 가져옴
        groups_info_list = config['group_info']
        rois_info_list = config['roi_info']
        
        expected_files = []
        
        # (수정) config의 'group_info' 리스트를 순회
        for group_dict in groups_info_list:
            group_name = group_dict['name']
            n_count = group_dict['n_per_group']
            s_count = group_dict['num_serial']
            
            n_range = range(1, n_count + 1)
            s_range = range(1, s_count + 1)
            
            # (수정) config의 'roi_info' 리스트를 순회
            for roi_dict in rois_info_list:
                roi_name = roi_dict['name']
                
                # n과 s의 조합 생성
                for n, s in itertools.product(n_range, s_range):
                    # 파일명 템플릿: "{name_of_group} {n}-{s} {name_of_roi}-total ion current.imzML"
                    file_name = f"{group_name} {n}-{s} {roi_name}-total ion count.imzML"
                    full_path = os.path.join(data_dir, file_name)
                    expected_files.append(full_path)
            
        return expected_files

    except KeyError as e:
        print(f"오류: config.yaml 파일에 필요한 키가 없습니다: {e}")
        return None
    except Exception as e:
        print(f"파일 목록 생성 중 오류 발생: {e}")
        return None

def validate_files(file_list):
    """
    (변경 없음)
    생성된 파일 목록이 실제 디스크에 모두 존재하는지 유효성 검사를 수행합니다.
    """
    print(f"\n--- 파일 유효성 검사 (총 {len(file_list)}개) ---")
    missing_files = []
    for f_path in file_list:
        if not os.path.exists(f_path):
            missing_files.append(f_path)
            
    if missing_files:
        print(f"오류: 다음 {len(missing_files)}개의 파일을 찾을 수 없습니다:")
        for missing in missing_files:
            print(f"  - {missing}")
        # .ibd 파일 누락 가능성도 안내
        print("\n.imzML 파일과 .ibd 파일이 모두 존재하는지 확인하세요.")
        return False
    
    print(f"성공: 모든 예상 파일 {len(file_list)}개를 data_dir에서 확인했습니다.")
    return True

def main(config_path='config/config.yaml'):
    """
    (수정됨) 메인 파이프라인 실행 함수
    config_path를 인자로 받아 해당 설정으로 파싱을 수행합니다.
    """
    print("--- 1. 설정 파일 로드 ---")
    # (수정) 하드코딩된 경로 대신 인자로 받은 config_path 사용
    config = parsing.load_yaml(config_path)
    if config is None:
        return # (Jupyter 호환)

    print("\n--- 2. m/z 빈(Bin) 정보 로드 ---")
    if 'binning_settings' not in config:
        print(f"오류: '{config_path}'에 'binning_settings' 섹션이 없습니다.")
        return # (Jupyter 호환)
        
    master_bins, bin_names = parsing.load_master_bins(config['binning_settings'])
    
    if master_bins is None or bin_names is None:
        print("m/z 빈 로드에 실패하여 처리를 중단합니다.")
        return # (Jupyter 호환)

    print("\n--- 3. 처리할 파일 목록 생성 ---")
    file_list = generate_expected_files(config)
    if file_list is None:
        print("파일 목록 생성에 실패하여 처리를 중단합니다.")
        return # (Jupyter 호환)

    # 4. 파일 유효성 검사
    if not validate_files(file_list):
        print("\n처리를 중단합니다. 파일 경로와 config.yaml 설정을 확인하세요.")
        return # (Jupyter 호환)
        
    # 5. 모든 파일 순차 처리
    print(f"\n--- 4. 총 {len(file_list)}개 파일 처리 시작 ---")
    output_dir = config.get('output_dir', 'output') # config에서 output_dir 읽기
    
    # 결과 저장 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    for i, imzml_file in enumerate(file_list):
        print(f"\n--- [ {i+1} / {len(file_list)} ] ---")
        # parsing.py의 함수를 호출하여 처리
        parsing.process_imzml_with_bins(imzml_file, master_bins, bin_names, output_dir)
        
    print(f"\n--- 모든 파일 처리가 완료되었습니다. ---")
    print(f"결과물은 '{output_dir}' 폴더에 저장되었습니다.")

if __name__ == '__main__':
    # 스크립트를 직접 실행할 경우 기본 config 파일을 사용
    main(config_path='config/config.yaml')