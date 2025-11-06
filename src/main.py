import os
import itertools
import yaml
# parsing.py에서 필요한 함수들을 import합니다.
try:
    import parsing
except ImportError:
    print("오류: 'parsing.py' 파일을 찾을 수 없습니다.")
    print("main.py와 parsing.py가 동일한 디렉토리에 있는지 확인하세요.")
    exit(1)

CONFIG_FILE = 'config.yaml'

def generate_expected_files(config):
    """
    config.yaml의 정보를 바탕으로 예상되는 모든 .imzML 파일의
    전체 경로 리스트를 생성합니다.
    """
    try:
        data_dir = config['data_dir']
        groups = config['group_info']['name_of_group']
        n_count = config['n_per_group']
        s_count = config['num_serial']
        rois = config['roi_info']['name_of_roi']
        
        # n과 s는 1부터 시작하는 숫자입니다.
        n_range = range(1, n_count + 1)
        s_range = range(1, s_count + 1)
        
        # 모든 조합 생성 (itertools.product 사용)
        all_combinations = itertools.product(groups, n_range, s_range, rois)
        
        expected_files = []
        # 파일명 템플릿: "{name_of_group} {n}-{s} {name_of_roi}-total ion current.imzML"
        # 참고: config.yaml의 예시(A_1-1...)와 파일명 템플릿이 다를 경우,
        # 이 템플릿 문자열을 실제 파일명 규칙에 맞게 수정해야 합니다.
        for (group, n, s, roi) in all_combinations:
            # 사용자님이 명시한 파일 이름 템플릿
            file_name = f"{group} {n}-{s} {roi}-total ion current.imzML"
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

def main():
    """
    메인 파이프라인 실행 함수
    """
    print("--- 1. 설정 파일 로드 ---")
    config = parsing.load_yaml(CONFIG_FILE)
    if config is None:
        exit(1)

    print("\n--- 2. m/z 빈(Bin) 정보 로드 ---")
    if 'binning_settings' not in config:
        print(f"오류: '{CONFIG_FILE}'에 'binning_settings' 섹션이 없습니다.")
        exit(1)
        
    master_bins, bin_names = parsing.load_master_bins(config['binning_settings'])
    
    if master_bins is None or bin_names is None:
        print("m/z 빈 로드에 실패하여 처리를 중단합니다.")
        exit(1)

    print("\n--- 3. 처리할 파일 목록 생성 ---")
    file_list = generate_expected_files(config)
    if file_list is None:
        print("파일 목록 생성에 실패하여 처리를 중단합니다.")
        exit(1)

    # 4. 파일 유효성 검사
    if not validate_files(file_list):
        print("\n처리를 중단합니다. 파일 경로와 config.yaml 설정을 확인하세요.")
        exit(1)
        
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
    main()
