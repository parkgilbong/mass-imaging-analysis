from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import pandas as pd
import os
import yaml  # PyYAML 라이브러리. 'pip install pyyaml' 필요

def load_yaml(config_path):
    """YAML 설정 파일을 로드합니다."""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print(f"오류: 설정 파일({config_path})을 찾을 수 없습니다.")
        return None
    except Exception as e:
        print(f"오류: 설정 파일 로드 중 문제 발생 - {e}")
        return None

def load_master_bins(binning_config):
    """
    config의 binning_settings를 기반으로 
    마스터 m/z 빈 리스트 [(min_mz, max_mz), ...]와 
    컬럼 이름 리스트 [bin_name, ...]를 반환합니다.
    """
    mode = binning_config.get('mode', 'file')
    master_bins = []
    bin_names = []
    
    print(f"m/z 빈 로드 모드: '{mode}'")

    try:
        if mode == 'file':
            file_path = binning_config['file_path']
            print(f"'{file_path}' 파일에서 m/z 빈 목록을 로드합니다.")
            
            # SCiLS에서 내보낸 'Mass ranges...' 파일 파싱
            # 주석('#')으로 시작하는 줄을 건너뛰고, 구분자는 세미콜론(;) 사용
            df = pd.read_csv(file_path, comment='#', delimiter=';')
            
            # 'm/z'와 'Interval Width (+/- Da)' 컬럼을 사용하여 min/max 범위 계산
            for _, row in df.iterrows():
                mz_center = row['m/z']
                mz_width = row['Interval Width (+/- Da)']
                min_mz = mz_center - mz_width
                max_mz = mz_center + mz_width
                
                # 'Name' 컬럼이 비어있으면 m/z 값으로 컬럼 이름 생성
                bin_name = row.get('Name')
                if pd.isna(bin_name) or bin_name.strip() == "":
                    bin_name = f"{mz_center:.4f}"
                
                master_bins.append((min_mz, max_mz))
                bin_names.append(bin_name)
                
        elif mode == 'direct':
            print("config.yaml의 'direct_bins'에서 m/z 빈 목록을 로드합니다.")
            for bin_info in binning_config['direct_bins']:
                mz_center = bin_info['mz']
                mz_width = bin_info['width']
                bin_name = bin_info.get('name', f"{mz_center:.4f}")
                
                min_mz = mz_center - mz_width
                max_mz = mz_center + mz_width
                
                master_bins.append((min_mz, max_mz))
                bin_names.append(bin_name)
        
        else:
            print(f"오류: 알 수 없는 binning_settings mode: '{mode}'")
            return None, None

        print(f"총 {len(master_bins)}개의 m/z 빈을 로드했습니다.")
        return master_bins, bin_names

    except FileNotFoundError:
        print(f"오류: m/z 빈 파일({binning_config.get('file_path')})을 찾을 수 없습니다.")
        return None, None
    except KeyError as e:
        print(f"오류: binning_settings에 필요한 키가 없습니다: {e}")
        return None, None
    except Exception as e:
        print(f"m/z 빈 로드 중 오류 발생: {e}")
        return None, None


def process_imzml_with_bins(imzml_filepath, master_bins, bin_names):
    """
    .imzML 파일을 파싱하여, 미리 정의된 'master_bins'를 기준으로
    intensity 매트릭스와 m/z 채널별 평균 intensity를 CSV 파일로 저장합니다.

    Args:
        imzml_filepath (str): 분석할 .imzML 파일의 경로.
        master_bins (list): (min_mz, max_mz) 튜플의 리스트.
        bin_names (list): CSV 컬럼 헤더로 사용할 m/z 빈 이름 리스트.
    """
    try:
        # 1. ImzMLParser 객체 생성
        p = ImzMLParser(imzml_filepath, parse_lib='lxml')
        print(f"\n파일 열기 성공: {imzml_filepath}")

        num_spectra = len(p.coordinates)
        num_bins = len(master_bins)
        
        if num_spectra == 0:
            print("오류: 파일에 스펙트럼 데이터가 없습니다.")
            return

        print(f"총 스펙트럼 수: {num_spectra}")
        print(f"적용할 m/z 빈(컬럼) 수: {num_bins}")

        # 2. 마스터 m/z 축(bins)에 맞게 각 스펙트럼의 intensity를 재구성합니다.
        all_aligned_intensities = []
        coordinates_list = []
        print("1단계: 각 스펙트럼을 마스터 m/z 빈에 정렬하는 중...")

        for i, (x, y, z) in enumerate(p.coordinates):
            mzs, intensities = p.getspectrum(i)
            
            # 마스터 m/z 축 길이에 맞는 0으로 채워진 배열 생성
            aligned_intensities = np.zeros(num_bins, dtype=np.float32)
            
            # 2-1. (핵심) 현재 스펙트럼의 각 m/z 피크를 마스터 빈과 대조
            for mz_val, intensity_val in zip(mzs, intensities):
                # 2-2. 이 m/z 값이 어느 마스터 빈에 속하는지 검색
                for bin_index, (min_mz, max_mz) in enumerate(master_bins):
                    if min_mz <= mz_val <= max_mz:
                        # 2-3. 해당 빈에 intensity를 합산
                        aligned_intensities[bin_index] += intensity_val
                        # 하나의 m/z 값은 하나의 빈에만 속한다고 가정
                        break 
            
            all_aligned_intensities.append(aligned_intensities)
            coordinates_list.append({'x': x, 'y': y})
            
            if (i + 1) % 1000 == 0 or (i + 1) == num_spectra:
                print(f"  ... {i + 1} / {num_spectra} 스펙트럼 처리 완료")

        print("데이터 정렬 완료.")

        # 3. pandas DataFrame 생성
        print("2단계: CSV 파일 생성 중...")
        df_coords = pd.DataFrame(coordinates_list)
        df_intensities = pd.DataFrame(all_aligned_intensities, columns=bin_names, dtype=np.float32)
        
        df_full = pd.concat([df_coords, df_intensities], axis=1)

        # 4. 전체 intensity 매트릭스를 CSV 파일로 저장
        base_filename = os.path.splitext(imzml_filepath)[0]
        output_intensities_csv = f"{base_filename}_binned_spectra.csv"
        df_full.to_csv(output_intensities_csv, index=False, float_format='%.4f')
        print(f"\nBinned Intensity 매트릭스 저장 완료: {output_intensities_csv}")
        print(f"  - 형태: {df_full.shape}")

        # 5. m/z 채널별 평균 intensity 계산 및 CSV 파일로 저장
        mean_intensities = df_intensities.mean().to_frame().T
        output_mean_csv = f"{base_filename}_mean_intensities.csv"
        mean_intensities.to_csv(output_mean_csv, index=False, float_format='%.4f')
        print(f"평균 Intensity 저장 완료: {output_mean_csv}")
        print(f"  - 형태: {mean_intensities.shape}")

    except FileNotFoundError:
        print(f"오류: 파일을 찾을 수 없습니다. {imzml_filepath}")
        print("(.imzML 파일과 .ibd 파일이 모두 존재하는지 확인하세요.)")
    except ImportError:
        print("오류: 'lxml' 라이브러리를 찾을 수 없습니다.")
        print("빠른 파싱을 위해 'pip install lxml' 명령어로 설치해주세요.")
    except Exception as e:
        print(f"데이터 파싱 또는 처리 중 오류 발생: {e}")


if __name__ == '__main__':
    # --- 실행 ---
    
    CONFIG_FILE = 'config.yaml'
    
    # 1. 설정 파일 로드
    config = load_yaml(CONFIG_FILE)
    if config is None:
        exit(1)

    # 2. 설정 파일에서 m/z 빈(bins) 정보 로드
    if 'binning_settings' not in config:
        print(f"오류: '{CONFIG_FILE}'에 'binning_settings' 섹션이 없습니다.")
        exit(1)
        
    master_bins, bin_names = load_master_bins(config['binning_settings'])
    
    if master_bins is None or bin_names is None:
        print("m/z 빈 로드에 실패하여 처리를 중단합니다.")
        exit(1)

    # 3. 분석할 .imzML 파일 경로 지정 (예시)
    # 실제 사용 시 이 부분을 config 파일에서 읽어오거나 인자로 받도록 수정할 수 있습니다.
    imzml_file = 'A_1-1_cortex-total_ion_count.imzML'

    if not os.path.exists(imzml_file):
        print(f"입력 파일 '{imzml_file}'을(를) 찾을 수 없습니다.")
        print("분석할 .imzML 파일의 경로를 'imzml_file' 변수에 올바르게 지정해주세요.")
    else:
        # 4. 처리 실행
        process_imzml_with_bins(imzml_file, master_bins, bin_names)
