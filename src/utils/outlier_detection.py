"""
Outlier detection utilities for mass imaging analysis.

Provides multiple methods for detecting and removing outliers from intensity data:
- Modified Z-score (MAD-based): Robust to non-normal distributions
- Standard Z-score: Assumes normal distribution
- IQR (Interquartile Range): Distribution-free method
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def detect_outliers_modified_z(data, cutoff=3.5):
    """
    Detect outliers using Modified Z-score (MAD-based).
    Robust to non-normal distributions.
    
    The modified Z-score uses the Median Absolute Deviation (MAD) instead of
    standard deviation, making it more robust to outliers in the data itself.
    
    Args:
        data: numpy array of values
        cutoff: threshold for modified z-score (default: 3.5)
    
    Returns:
        Boolean mask where True = outlier
    """
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    
    if mad == 0:
        # All values are identical or very close, no outliers
        return np.zeros(len(data), dtype=bool)
    
    # 0.6745 is the constant to convert MAD to standard deviation equivalent
    # for normally distributed data
    modified_z_scores = 0.6745 * (data - median) / mad
    return np.abs(modified_z_scores) > cutoff


def detect_outliers_z_score(data, cutoff=3.0):
    """
    Detect outliers using standard Z-score.
    Assumes normal distribution.
    
    Args:
        data: numpy array of values
        cutoff: threshold for z-score (default: 3.0)
    
    Returns:
        Boolean mask where True = outlier
    """
    mean = np.mean(data)
    std = np.std(data)
    
    if std == 0:
        # All values are identical, no outliers
        return np.zeros(len(data), dtype=bool)
    
    z_scores = (data - mean) / std
    return np.abs(z_scores) > cutoff


def detect_outliers_iqr(data, cutoff=1.5):
    """
    Detect outliers using Interquartile Range (IQR).
    Very robust, distribution-free method.
    
    Outliers are defined as values below Q1 - cutoff*IQR or above Q3 + cutoff*IQR.
    
    Args:
        data: numpy array of values
        cutoff: multiplier for IQR (default: 1.5, standard Tukey's fences)
    
    Returns:
        Boolean mask where True = outlier
    """
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    
    if iqr == 0:
        # All values in the interquartile range are identical
        return np.zeros(len(data), dtype=bool)
    
    lower_bound = q1 - cutoff * iqr
    upper_bound = q3 + cutoff * iqr
    
    return (data < lower_bound) | (data > upper_bound)


def remove_outliers(data, method='modified_z', cutoff=3.5, min_remaining=3):
    """
    Remove outliers from data array using specified method.
    
    Args:
        data: numpy array of values
        method: 'modified_z', 'z_score', or 'iqr'
        cutoff: threshold value (interpretation depends on method)
        min_remaining: minimum points to keep (safety check)
    
    Returns:
        tuple: (filtered_data, n_removed, outlier_mask)
            - filtered_data: array with outliers removed
            - n_removed: number of outliers removed
            - outlier_mask: boolean array indicating outliers
    """
    if len(data) < min_remaining:
        logger.warning(
            f"Data size ({len(data)}) < min_remaining ({min_remaining}), "
            f"skipping outlier removal"
        )
        return data, 0, np.zeros(len(data), dtype=bool)
    
    # Detect outliers based on method
    if method == 'modified_z':
        outlier_mask = detect_outliers_modified_z(data, cutoff)
    elif method == 'z_score':
        outlier_mask = detect_outliers_z_score(data, cutoff)
    elif method == 'iqr':
        outlier_mask = detect_outliers_iqr(data, cutoff)
    else:
        logger.error(f"Unknown outlier detection method: '{method}', using original data")
        return data, 0, np.zeros(len(data), dtype=bool)
    
    n_outliers = np.sum(outlier_mask)
    n_remaining = len(data) - n_outliers
    
    # Safety check: ensure minimum remaining points
    if n_remaining < min_remaining:
        logger.warning(
            f"Outlier removal would leave {n_remaining} points (< {min_remaining}), "
            f"keeping all {len(data)} points"
        )
        return data, 0, np.zeros(len(data), dtype=bool)
    
    filtered_data = data[~outlier_mask]
    
    if n_outliers > 0:
        logger.debug(
            f"Removed {n_outliers}/{len(data)} outliers "
            f"({100*n_outliers/len(data):.1f}%)"
        )
    
    return filtered_data, n_outliers, outlier_mask
