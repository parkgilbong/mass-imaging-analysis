from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import pandas as pd
import os

def process_imzml_to_csv(imzml_filepath, mz_precision=1):
    """
    .imzML 파일을 파싱하여 모든 스펙트럼의 intensity 매트릭스와
    m/z 채널별 평균 intensity를 각각 CSV 파일로 저장합니다.

    Args:
        imzml_filepath (str): 분석할 .imzML 파일의 경로.
        mz_precision (int): m/z 값을 그룹화할 소수점 자릿수.
    """
    try:
        # 1. ImzMLParser 객체 생성 (lxml 라이브러리로 빠른 파싱)
        p = ImzMLParser(imzml_filepath, parse_lib='lxml')
        print(f"파일 열기 성공: {imzml_filepath}")

        num_spectra = len(p.coordinates)
        if num_spectra == 0:
            print("오류: 파일에 스펙트럼 데이터가 없습니다.")
            return

        print(f"총 스펙트럼 수: {num_spectra}")

        # 2. 모든 스펙트럼을 스캔하여 고유한 m/z 값의 "마스터 리스트"를 생성합니다.
        print(f"1단계: m/z 값을 소수점 {mz_precision}자리까지 반올림하여 마스터 축 생성 중...")
        master_mzs_set = set()
        for i, _ in enumerate(p.coordinates):
            mzs, _ = p.getspectrum(i)
            master_mzs_set.update(np.round(mzs, mz_precision))
        
        master_mzs = sorted(list(master_mzs_set))
        num_master_mzs = len(master_mzs)
        mz_to_index_map = {mz: i for i, mz in enumerate(master_mzs)}
        print(f"  - 고유 m/z 값 {num_master_mzs}개를 기준으로 데이터 테이블을 생성합니다.")

        # 3. 마스터 m/z 축에 맞게 각 스펙트럼의 intensity를 재구성합니다.
        all_aligned_intensities = []
        coordinates_list = []
        print("2단계: 각 스펙트럼을 마스터 m/z 축에 정렬하는 중...")
        for i, (x, y, z) in enumerate(p.coordinates):
            mzs, intensities = p.getspectrum(i)
            
            # 마스터 m/z 축 길이에 맞는 0으로 채워진 배열 생성
            aligned_intensities = np.zeros(num_master_mzs, dtype=np.float32)
            
            # 현재 스펙트럼의 m/z 값을 마스터 축과 동일한 정밀도로 반올림
            rounded_mzs = np.round(mzs, mz_precision)

            # 현재 스펙트럼의 m/z 값을 마스터 축의 인덱스와 매핑하여 intensity 값을 할당
            for mz, intensity in zip(rounded_mzs, intensities):
                if mz in mz_to_index_map:
                    index = mz_to_index_map[mz]
                    aligned_intensities[index] += intensity # 동일한 bin에 속하는 intensity는 합산
            
            all_aligned_intensities.append(aligned_intensities)
            coordinates_list.append({'x': x, 'y': y})
        print("데이터 읽기 완료.")

        # 4. pandas DataFrame 생성
        # 헤더(컬럼명)를 소수점 둘째 자리까지 포맷팅
        formatted_headers = [f"{mz:.2f}" for mz in master_mzs]
        df_coords = pd.DataFrame(coordinates_list)
        df_intensities = pd.DataFrame(all_aligned_intensities, columns=formatted_headers)
        
        # 좌표와 intensity 데이터를 합칩니다.
        df_full = pd.concat([df_coords, df_intensities], axis=1)

        # 5. 전체 intensity 매트릭스를 CSV 파일로 저장
        base_filename = os.path.splitext(imzml_filepath)[0]
        output_intensities_csv = f"{base_filename}_full_spectra.csv"
        df_full.to_csv(output_intensities_csv, index=False, float_format='%.4f')
        print(f"\nIntensity 매트릭스 저장 완료: {output_intensities_csv}")
        print(f"  - 형태: {df_full.shape}")

        # 6. m/z 채널별 평균 intensity 계산 및 CSV 파일로 저장
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
    # --- 사용 예시 ---

    # 분석할 .imzML 파일 경로
    # .ibd 파일도 동일한 경로에 있어야 합니다.
    imzml_file = 'A_1-1_cortex-total_ion_count.imzML'

    # 파일이 존재하는지 확인
    if not os.path.exists(imzml_file):
        print(f"입력 파일 '{imzml_file}'을(를) 찾을 수 없습니다.")
        print("분석할 .imzML 파일의 경로를 'imzml_file' 변수에 올바르게 지정해주세요.")
        # 예시용 더미 파일 생성 (실제 분석 시에는 이 부분을 제거하거나 주석 처리)
        print("테스트를 위해 더미 파일을 생성합니다.")
        from pyimzml.ImzMLWriter import ImzMLWriter
        mzs = np.linspace(100, 1000, 7)
        coords = [(x, y, 1) for x in range(32) for y in range(32)] # 1024 spectra
        coords.append((32,1,1)) # 1025 spectra
        with ImzMLWriter(imzml_file, mode="w") as writer:
            for x, y, z in coords:
                intensities = np.random.rand(len(mzs)) * 1000
                writer.addSpectrum(mzs, intensities, (x, y, z))
        print(f"더미 파일 '{imzml_file}' 생성 완료.")

    # 함수 호출하여 처리 실행
    process_imzml_to_csv(imzml_file, mz_precision=2)