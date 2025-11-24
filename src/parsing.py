import warnings
# pyimzml에서 발생하는 특정 UserWarning을 무시합니다.
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')

from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import pandas as pd
import os
import yaml

# (NEW) 로깅 유틸리티 임포트
try:
    from utils.logging_utils import get_logger
except ImportError:
    try:
        from .utils.logging_utils import get_logger
    except ImportError:
        import logging
        get_logger = logging.getLogger

# 로거 초기화
logger = get_logger(__name__)

def load_yaml(config_path):
    """YAML 설정 파일을 로드합니다."""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"설정 파일({config_path})을 찾을 수 없습니다.")
        return None
    except Exception as e:
        logger.error(f"설정 파일 로드 중 문제 발생: {e}")
        return None

def load_master_bins(binning_config):
    """
    config의 binning_settings를 기반으로 m/z 빈 리스트를 반환합니다.
    """
    mode = binning_config.get('mode', 'file')
    master_bins = []
    bin_names = []
    
    logger.info(f"m/z bin load 모드: '{mode}'")

    try:
        if mode == 'file':
            file_path = binning_config['file_path']
            logger.info(f"'{file_path}' 파일에서 m/z bin 목록을 로드합니다.")
            
            if not os.path.exists(file_path):
                 logger.error(f"파일을 찾을 수 없습니다: {file_path}")
                 return None, None

            df = pd.read_csv(file_path, comment='#', delimiter=';')
            
            for _, row in df.iterrows():
                mz_center = row['m/z']
                mz_width = row['Interval Width (+/- Da)']
                min_mz = mz_center - mz_width
                max_mz = mz_center + mz_width
                
                bin_name = row.get('Name')
                if pd.isna(bin_name) or (isinstance(bin_name, str) and bin_name.strip() == ""):
                    bin_name = f"{mz_center:.4f}"
                
                master_bins.append((min_mz, max_mz))
                bin_names.append(bin_name)
                
        elif mode == 'direct':
            logger.info("config.yaml의 'direct_bins'에서 m/z bin 목록을 로드합니다.")
            for bin_info in binning_config['direct_bins']:
                mz_center = bin_info['mz']
                mz_width = bin_info['width']
                bin_name = bin_info.get('name', f"{mz_center:.4f}")
                
                min_mz = mz_center - mz_width
                max_mz = mz_center + mz_width
                
                master_bins.append((min_mz, max_mz))
                bin_names.append(bin_name)
        else:
            logger.error(f"알 수 없는 binning_settings mode: '{mode}'")
            return None, None

        logger.info(f"총 {len(master_bins)}개의 m/z bin을 로드했습니다.")
        return master_bins, bin_names

    except Exception as e:
        logger.error(f"m/z 빈 로드 중 오류 발생: {e}", exc_info=True)
        return None, None


def process_imzml_with_bins(imzml_filepath, master_bins, bin_names, output_dir="."):
    """
    .imzML 파일을 파싱하여 CSV 파일로 저장합니다.
    """
    try:
        # 1. ImzMLParser 객체 생성
        # (parse_lib='xml' 사용 유지 - 호환성 위함)
        p = ImzMLParser(imzml_filepath, parse_lib='xml')
        logger.info(f"파일 열기 성공: {os.path.basename(imzml_filepath)}")

        num_spectra = len(p.coordinates)
        num_bins = len(master_bins)
        
        if num_spectra == 0:
            logger.warning(f"파일에 스펙트럼 데이터가 없습니다: {imzml_filepath}")
            return

        # logger.info(f"총 스펙트럼 수: {num_spectra}, 적용할 빈 수: {num_bins}")

        # 2. 마스터 m/z 축(bins)에 맞게 각 스펙트럼의 intensity를 재구성
        all_aligned_intensities = []
        coordinates_list = []
        
        # 로깅 양 조절을 위해 10% 단위나 1000개 단위로 로그 출력
        log_interval = max(1, num_spectra // 5) 

        for i, (x, y, z) in enumerate(p.coordinates):
            mzs, intensities = p.getspectrum(i)
            
            aligned_intensities = np.zeros(num_bins, dtype=np.float32)
            
            for bin_index, (min_mz, max_mz) in enumerate(master_bins):
                # 현재 스펙트럼에서 해당 bin 범위에 있는 intensity 합산
                # (속도 최적화를 위해 numpy boolean indexing 대신 루프나 searchsorted를 쓸 수도 있지만
                #  현재는 가독성을 위해 유지)
                mask = (mzs >= min_mz) & (mzs <= max_mz)
                if mask.any():
                    aligned_intensities[bin_index] = np.sum(intensities[mask])
            
            all_aligned_intensities.append(aligned_intensities)
            coordinates_list.append({'x': x, 'y': y})
            
            if (i + 1) % log_interval == 0:
                logger.info(f"  ... {i + 1}/{num_spectra} 스펙트럼 처리 중")

        # 3. pandas DataFrame 생성
        df_coords = pd.DataFrame(coordinates_list)
        df_intensities = pd.DataFrame(all_aligned_intensities, columns=bin_names, dtype=np.float32)
        df_full = pd.concat([df_coords, df_intensities], axis=1)

        # 4. 저장
        base_filename = os.path.splitext(os.path.basename(imzml_filepath))[0]
        output_intensities_csv = os.path.join(output_dir, f"{base_filename}_binned_spectra.csv")
        df_full.to_csv(output_intensities_csv, index=False, float_format='%.4f')
        # logger.info(f"Binned Intensity 저장 완료: {output_intensities_csv}")

        # 5. 평균 저장
        mean_intensities = df_intensities.mean().to_frame().T
        output_mean_csv = os.path.join(output_dir, f"{base_filename}_mean_intensities.csv")
        mean_intensities.to_csv(output_mean_csv, index=False, float_format='%.4f')
        logger.info(f"처리 완료 및 저장됨: {output_mean_csv}")

    except FileNotFoundError:
        logger.error(f"파일을 찾을 수 없습니다: {imzml_filepath}")
    except ImportError:
        logger.error("'lxml' 라이브러리 오류.")
    except Exception as e:
        logger.error(f"데이터 처리 중 오류 발생 ({os.path.basename(imzml_filepath)}): {e}", exc_info=True)