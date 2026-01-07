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
# Logger will be initialized in main() with optional log file

def _parse_individuals_info(group_dict):
    """
    Parse per-individual serial section information from group configuration
    
    Args:
        group_dict: Group configuration dictionary
    
    Returns:
        List[Tuple[int, List[int]]]: List in format [(n, [s1, s2, ...]), ...]
    
    Raises:
        ValueError: If config format is invalid or validation fails
    """
    group_name = group_dict.get('name', 'unknown')
    
    n_per_group = group_dict.get('n_per_group')
    num_serial = group_dict.get('num_serial')
    
    # Determine individual IDs
    individual_ids = group_dict.get('individual_ids')
    if individual_ids:
        if not isinstance(individual_ids, list):
             raise ValueError(f"Group '{group_name}': individual_ids must be a list of integers.")
        
        # Infer n_per_group if missing
        if n_per_group is None:
            n_per_group = len(individual_ids)
        elif len(individual_ids) != n_per_group:
            raise ValueError(
                f"Group '{group_name}': individual_ids length ({len(individual_ids)}) "
                f"does not match n_per_group ({n_per_group}). "
                f"individual_ids: {individual_ids}"
            )
    else:
        if n_per_group is None:
             raise ValueError(f"Group '{group_name}': 'n_per_group' is required if 'individual_ids' is missing")
        individual_ids = list(range(1, n_per_group + 1))

    # Check for explicit serial_ids (overrides num_serial)
    serial_ids = group_dict.get('serial_ids')
    if serial_ids:
        if len(serial_ids) != n_per_group:
            raise ValueError(
                f"Group '{group_name}': serial_ids length ({len(serial_ids)}) "
                f"does not match n_per_group ({n_per_group}). "
                f"serial_ids: {serial_ids}"
            )
        # Validate structure (list of lists of ints)
        if not all(isinstance(s_list, list) and all(isinstance(s, int) for s in s_list) for s_list in serial_ids):
             raise ValueError(f"Group '{group_name}': serial_ids must be a list of lists of integers.")
        
        return list(zip(individual_ids, serial_ids))

    if num_serial is None:
        raise ValueError(f"Group '{group_name}': 'num_serial' (or 'serial_ids') is required")

    # Fallback: Generate serial IDs from num_serial
    # Case 1: num_serial is integer (traditional - all individuals same)
    if isinstance(num_serial, int):
        return [(ind_id, list(range(1, num_serial + 1))) for ind_id in individual_ids]
    
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
        
        return [(ind_id, list(range(1, count + 1))) for ind_id, count in zip(individual_ids, num_serial)]
    
    else:
        raise ValueError(
            f"Group '{group_name}': num_serial must be an integer or list of integers. "
            f"Current type: {type(num_serial)}"
        )

def main(config_path='config/config.yaml'):
    """Data aggregation execution function"""
    logger.info("========== Step 2: Starting data aggregation ==========")
    
    config = parsing.load_yaml(config_path)
    if config is None:
        return

    try:
        output_dir = config['output_dir']
        groups_info_list = config['group_info']
        rois_info_list = config['roi_info']
    except KeyError as e:
        logger.error(f"config.yaml key error: {e}")
        return
        
    logger.info(f"Aggregating data from '{output_dir}' folder.")

    all_aggregated_data = []
    m_z_columns = None

    # Calculate expected number of combinations
    # Note: This calculation is approximate now due to flexible serial IDs
    # total_combinations = len(groups_info_list) * sum(g['n_per_group'] for g in groups_info_list) * len(rois_info_list)
    # logger.info(f"Attempting to aggregate {total_combinations} combinations (group*n*roi)")

    for group_dict in groups_info_list:
        group_name = group_dict['name']
        
        # Parse per-individual serial section info (includes validation)
        try:
            individuals_info = _parse_individuals_info(group_dict)
        except ValueError as e:
            logger.error(f"Config validation error: {e}")
            return
        
        for roi_dict in rois_info_list:
            roi_name = roi_dict['name']
            
            # Process only actual serial section count per individual
            for n, serial_id_list in individuals_info:
                serial_data_to_average = []
                
                for s in serial_id_list:
                    base_imzml_name = f"{group_name} {n}-{s} {roi_name}-total ion count"
                    mean_csv_path = os.path.join(output_dir, f"{base_imzml_name}_mean_intensities.csv")
                    
                    if os.path.exists(mean_csv_path):
                        try:
                            df_mean = pd.read_csv(mean_csv_path)
                            serial_data_to_average.append(df_mean)
                            
                            if m_z_columns is None:
                                m_z_columns = df_mean.columns.tolist()
                                # logger.info(f"Detected m/z columns: {len(m_z_columns)} columns")
                        except Exception as e:
                            logger.error(f"Failed to read file ({mean_csv_path}): {e}")
                    else:
                        logger.warning(f"File missing: {mean_csv_path}")

                if not serial_data_to_average:
                    logger.warning(f"No data - skipping: {group_name} n={n} {roi_name}")
                    continue
                    
                df_concat = pd.concat(serial_data_to_average)
                df_averaged = df_concat.mean(axis=0)
                
                agg_row_data = df_averaged.to_dict()
                agg_row_data['group'] = group_name
                agg_row_data['n'] = n
                agg_row_data['roi'] = roi_name
                
                all_aggregated_data.append(agg_row_data)
                logger.info(f"Aggregation complete: {group_name} n={n} {roi_name} ({len(serial_data_to_average)} files)")

    if not all_aggregated_data:
        logger.error("No aggregated data. Please verify that Step 1 was completed successfully.")
        return
        
    df_final = pd.DataFrame(all_aggregated_data)
    
    if m_z_columns:
        identifier_cols = ['group', 'n', 'roi']
        final_cols = identifier_cols + m_z_columns
        final_cols = [col for col in final_cols if col in df_final.columns] 
        df_final = df_final[final_cols]

    output_csv_path = os.path.join(output_dir, "aggregated_mean_intensities.csv")
    df_final.to_csv(output_csv_path, index=False, float_format='%.4f')
    
    logger.info(f"Final aggregated file saved: {output_csv_path}")
    logger.info("========== Step 2 complete ==========")

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Step 2: Aggregate mean intensities across technical replicates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Run with default config
  python src/aggregate.py
  
  # Run with custom config
  python src/aggregate.py --config config/custom_config.yaml
  
  # Run with Snakemake (custom log file)
  python src/aggregate.py --config config/config.yaml --log-file output/logs/aggregate_data.log
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