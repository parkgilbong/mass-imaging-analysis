import warnings
# Ignore specific UserWarnings from pyimzml
warnings.filterwarnings('ignore', category=UserWarning, module='pyimzml')

from pyimzml.ImzMLParser import ImzMLParser
import numpy as np
import pandas as pd
import os
import yaml

# Import logging utility
try:
    from utils.logging_utils import get_logger
    from utils.outlier_detection import remove_outliers
except ImportError:
    try:
        from .utils.logging_utils import get_logger
        from .utils.outlier_detection import remove_outliers
    except ImportError:
        import logging
        get_logger = logging.getLogger
        # Fallback: define dummy remove_outliers if import fails
        def remove_outliers(data, **kwargs):
            return data, 0, None

# Initialize logger
logger = get_logger(__name__)

def load_yaml(config_path):
    """Load YAML configuration file."""
    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        return None
    except Exception as e:
        logger.error(f"Error loading configuration file: {e}")
        return None

def load_master_bins(binning_config):
    """
    Return m/z bin list based on binning_settings in config.
    """
    mode = binning_config.get('mode', 'file')
    master_bins = []
    bin_names = []
    
    logger.info(f"m/z bin loading mode: '{mode}'")

    try:
        if mode == 'file':
            file_path = binning_config['file_path']
            logger.info(f"Loading m/z bin list from file: '{file_path}'")
            
            if not os.path.exists(file_path):
                 logger.error(f"File not found: {file_path}")
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
            logger.info("Loading m/z bin list from 'direct_bins' in config.yaml")
            for bin_info in binning_config['direct_bins']:
                mz_center = bin_info['mz']
                mz_width = bin_info['width']
                bin_name = bin_info.get('name', f"{mz_center:.4f}")
                
                min_mz = mz_center - mz_width
                max_mz = mz_center + mz_width
                
                master_bins.append((min_mz, max_mz))
                bin_names.append(bin_name)
        else:
            logger.error(f"Unknown binning_settings mode: '{mode}'")
            return None, None

        logger.info(f"Loaded {len(master_bins)} m/z bins")
        return master_bins, bin_names

    except Exception as e:
        logger.error(f"Error loading m/z bins: {e}", exc_info=True)
        return None, None


def process_imzml_with_bins(imzml_filepath, master_bins, bin_names, output_dir=".", outlier_config=None):
    """
    Parse .imzML file and save as CSV file.
    
    Args:
        imzml_filepath: Path to .imzML file
        master_bins: List of (min_mz, max_mz) tuples
        bin_names: List of bin names
        output_dir: Output directory path
        outlier_config: Dict with outlier detection settings (optional)
                       Keys: 'enabled', 'method', 'cutoff', 'min_remaining'
    """
    try:
        # 1. Create ImzMLParser object
        # (Keep using parse_lib='xml' for compatibility)
        p = ImzMLParser(imzml_filepath, parse_lib='xml')
        logger.info(f"Successfully opened file: {os.path.basename(imzml_filepath)}")

        num_spectra = len(p.coordinates)
        num_bins = len(master_bins)
        
        if num_spectra == 0:
            logger.warning(f"No spectrum data in file: {imzml_filepath}")
            return

        # logger.info(f"Total spectra: {num_spectra}, bins to apply: {num_bins}")

        # 2. Reconstruct intensity for each spectrum according to master m/z bins
        all_aligned_intensities = []
        coordinates_list = []
        
        # Log at 20% intervals to reduce logging volume
        log_interval = max(1, num_spectra // 5) 

        for i, (x, y, z) in enumerate(p.coordinates):
            mzs, intensities = p.getspectrum(i)
            
            aligned_intensities = np.zeros(num_bins, dtype=np.float32)
            
            for bin_index, (min_mz, max_mz) in enumerate(master_bins):
                # Sum intensities within the bin range from current spectrum
                # (Could use searchsorted for speed optimization instead of boolean indexing,
                #  but keeping current approach for readability)
                mask = (mzs >= min_mz) & (mzs <= max_mz)
                if mask.any():
                    aligned_intensities[bin_index] = np.sum(intensities[mask])
            
            all_aligned_intensities.append(aligned_intensities)
            coordinates_list.append({'x': x, 'y': y})
            
            if (i + 1) % log_interval == 0:
                logger.info(f"  ... Processing {i + 1}/{num_spectra} spectra")

        # 3. Create pandas DataFrame
        df_coords = pd.DataFrame(coordinates_list)
        df_intensities = pd.DataFrame(all_aligned_intensities, columns=bin_names, dtype=np.float32)
        df_full = pd.concat([df_coords, df_intensities], axis=1)

        # 4. Save binned spectra
        base_filename = os.path.splitext(os.path.basename(imzml_filepath))[0]
        output_intensities_csv = os.path.join(output_dir, f"{base_filename}_binned_spectra.csv")
        df_full.to_csv(output_intensities_csv, index=False, float_format='%.4f')
        # logger.info(f"Binned intensity saved: {output_intensities_csv}")

        # 5. Calculate mean intensities with optional outlier removal
        if outlier_config and outlier_config.get('enabled', False):
            mean_intensities = calculate_mean_with_outlier_removal(
                df_intensities, 
                outlier_config,
                logger
            )
        else:
            mean_intensities = df_intensities.mean().to_frame().T
        
        output_mean_csv = os.path.join(output_dir, f"{base_filename}_mean_intensities.csv")
        mean_intensities.to_csv(output_mean_csv, index=False, float_format='%.4f')
        logger.info(f"Processing complete and saved: {output_mean_csv}")

    except FileNotFoundError:
        logger.error(f"File not found: {imzml_filepath}")
    except ImportError:
        logger.error("'lxml' library error.")
    except Exception as e:
        logger.error(f"Error during data processing ({os.path.basename(imzml_filepath)}): {e}", exc_info=True)


def calculate_mean_with_outlier_removal(df_intensities, outlier_config, logger):
    """
    Calculate mean intensities with outlier detection and removal.
    
    Processes each bin (column) separately to detect and remove outliers
    before computing the mean intensity for that bin.
    
    Args:
        df_intensities: DataFrame with intensity values (rows=spectra, cols=bins)
        outlier_config: Dict with 'method', 'cutoff', 'min_remaining'
        logger: Logger instance
    
    Returns:
        DataFrame with mean intensities (1 row, columns=bins)
    """
    method = outlier_config.get('method', 'modified_z')
    cutoff = outlier_config.get('cutoff', 3.5)
    min_remaining = outlier_config.get('min_remaining', 3)
    
    mean_values = []
    total_outliers = 0
    total_points = 0
    
    # Process each bin (column) separately
    for col in df_intensities.columns:
        bin_data = df_intensities[col].values
        total_points += len(bin_data)
        
        # Remove outliers from this bin
        filtered_data, n_removed, _ = remove_outliers(
            bin_data, 
            method=method, 
            cutoff=cutoff, 
            min_remaining=min_remaining
        )
        
        total_outliers += n_removed
        mean_values.append(np.mean(filtered_data))
    
    if total_outliers > 0:
        logger.info(
            f"Outlier detection: removed {total_outliers}/{total_points} points "
            f"({100*total_outliers/total_points:.2f}%) using {method} method (cutoff={cutoff})"
        )
    else:
        logger.info(f"Outlier detection: no outliers detected using {method} method (cutoff={cutoff})")
    
    # Create DataFrame with same structure as original mean calculation
    return pd.DataFrame([mean_values], columns=df_intensities.columns)