import os
import itertools
import yaml

# Import logging utility
try:
    import parsing
    from utils.logging_utils import get_logger
except ImportError:
    try:
        from . import parsing
        from .utils.logging_utils import get_logger
    except ImportError:
        print("Error: Module import failed.")
        pass 

CONFIG_FILE = 'config/config.yaml'
# Logger will be initialized in main() with optional log file

def _parse_individuals_info(group_dict):
    """
    Parse per-individual serial section information from group configuration
    
    Args:
        group_dict: Group configuration dictionary
    
    Returns:
        List[Tuple[int, int]]: List in format [(n, num_serial), ...]
    
    Raises:
        ValueError: If config format is invalid or validation fails
    """
    group_name = group_dict.get('name', 'unknown')
    
    if 'n_per_group' not in group_dict:
        raise ValueError(f"Group '{group_name}': 'n_per_group' is required")
    
    n_per_group = group_dict['n_per_group']
    num_serial = group_dict.get('num_serial')
    
    if num_serial is None:
        raise ValueError(f"Group '{group_name}': 'num_serial' is required")
    
    # Case 1: num_serial is integer (traditional - all individuals same)
    if isinstance(num_serial, int):
        return [(n, num_serial) for n in range(1, n_per_group + 1)]
    
    # Case 2: num_serial is list (new - per-individual)
    elif isinstance(num_serial, list):
        # Validation: list length must match n_per_group
        if len(num_serial) != n_per_group:
            raise ValueError(
                f"Group '{group_name}': num_serial list length ({len(num_serial)}) "
                f"does not match n_per_group ({n_per_group}). "
                f"num_serial: {num_serial}"
            )
        
        # Verify all elements are positive integers
        if not all(isinstance(s, int) and s > 0 for s in num_serial):
            raise ValueError(
                f"Group '{group_name}': All values in num_serial list must be positive integers. "
                f"num_serial: {num_serial}"
            )
        
        return [(n+1, serial_count) for n, serial_count in enumerate(num_serial)]
    
    else:
        raise ValueError(
            f"Group '{group_name}': num_serial must be an integer or list of integers. "
            f"Current type: {type(num_serial)}"
        )

def generate_expected_files(config):
    """
    Generate a list of all expected .imzML file paths based on the flexible
    'dictionary list' structure in config.yaml.
    Supports different numbers of serial sections per individual.
    """
    try:
        data_dir = config['data_dir']
        groups_info_list = config['group_info']
        rois_info_list = config['roi_info']
        
        expected_files = []
        
        for group_dict in groups_info_list:
            group_name = group_dict['name']
            
            # Parse per-individual serial section info (includes validation)
            individuals_info = _parse_individuals_info(group_dict)
            
            for roi_dict in rois_info_list:
                roi_name = roi_dict['name']
                
                # Generate files only for actual serial section count per individual
                for n, num_serial in individuals_info:
                    for s in range(1, num_serial + 1):
                        file_name = f"{group_name} {n}-{s} {roi_name}-total ion count.imzML"
                        full_path = os.path.join(data_dir, file_name)
                        expected_files.append(full_path)
            
        return expected_files

    except KeyError as e:
        logger.error(f"config.yaml key error: {e}")
        return None
    except ValueError as e:
        logger.error(f"config.yaml validation error: {e}")
        return None
    except Exception as e:
        logger.error(f"Error generating file list: {e}")
        return None

def validate_files(file_list):
    logger.info(f"Starting file validation (expecting {len(file_list)} files total)")
    missing_files = []
    for f_path in file_list:
        if not os.path.exists(f_path):
            missing_files.append(f_path)
            
    if missing_files:
        logger.error(f"Cannot find the following {len(missing_files)} files:")
        for missing in missing_files:
            logger.error(f"  - {missing}")
        logger.error("Please verify that both .imzML and .ibd files exist.")
        return False
    
    logger.info("All expected files exist.")
    return True

def main(config_path='config/config.yaml'):
    """Main pipeline execution function"""
    logger.info("========== Step 1: Starting data parsing ==========")
    
    config = parsing.load_yaml(config_path)
    if config is None:
        return

    if 'binning_settings' not in config:
        logger.error("Config file is missing 'binning_settings'.")
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
    
    logger.info(f"Starting to process {len(file_list)} files total.")
    
    # Load outlier detection config
    outlier_config = config.get('outlier_detection', {'enabled': False})
    
    for i, imzml_file in enumerate(file_list):
        logger.info(f"--- Processing file [{i+1}/{len(file_list)}] ---")
        parsing.process_imzml_with_bins(
            imzml_file, 
            master_bins, 
            bin_names, 
            output_dir,
            outlier_config=outlier_config
        )
        
    logger.info("========== Step 1 complete: All files processed ==========")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Step 1: Parse imzML files and extract m/z bin intensities',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Run with default config
  python src/main.py
  
  # Run with custom config
  python src/main.py --config config/custom_config.yaml
  
  # Run with Snakemake (custom log file)
  python src/main.py --config config/config.yaml --log-file output/logs/parse_data.log
        '''
    )
    parser.add_argument(
        '--config',
        type=str,
        default='config/config.yaml',
        help='Path to config YAML file (default: config/config.yaml)'
    )
    parser.add_argument(
        '--log-file',
        type=str,
        default=None,
        help='Path to log file (for Snakemake integration)'
    )
    
    args = parser.parse_args()
    
    # Initialize logger with custom log file if provided
    global logger
    logger = get_logger(__name__, log_file=args.log_file)
    
    main(config_path=args.config)