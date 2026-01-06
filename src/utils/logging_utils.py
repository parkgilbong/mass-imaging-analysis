import logging
import os
import sys
from datetime import datetime

def get_logger(name, log_dir='logs', log_file=None):
    """
    Create and return a logger that logs to both console and file.
    
    Args:
        name (str): Logger name (usually __name__)
        log_dir (str): Directory to save log files (used if log_file is None)
        log_file (str): Specific log file path (overrides log_dir if provided)
    
    Returns:
        logging.Logger: Configured logger object
    """
    # 로거 생성
    logger = logging.getLogger(name)
    
    # 이미 핸들러가 설정되어 있다면(중복 방지) 그대로 반환
    if logger.handlers:
        return logger
    
    logger.setLevel(logging.INFO)
    
    # 포맷터 설정
    # 예: [2024-01-01 12:00:00] INFO [main.py:45] 처리 시작...
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)d] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # 1. 콘솔 핸들러 (StreamHandler)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.INFO)
    logger.addHandler(console_handler)
    
    # 2. File handler (FileHandler)
    try:
        if log_file:
            # Use specified log file (for Snakemake integration)
            log_filepath = log_file
            # Create directory if it doesn't exist
            log_file_dir = os.path.dirname(log_filepath)
            if log_file_dir and not os.path.exists(log_file_dir):
                os.makedirs(log_file_dir)
        else:
            # Use default date-based log file
            if not os.path.exists(log_dir):
                os.makedirs(log_dir)
            today = datetime.now().strftime('%Y-%m-%d')
            log_filepath = os.path.join(log_dir, f"{today}_analysis.log")
        
        file_handler = logging.FileHandler(log_filepath, encoding='utf-8')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.INFO)  # Change to DEBUG if needed
        logger.addHandler(file_handler)
        
    except Exception as e:
        print(f"[Warning] 로그 파일 생성 실패: {e}")
    
    # 전파 방지 (Jupyter Notebook 중복 출력 방지)
    logger.propagate = False
    
    return logger