import logging
import os
import sys
from datetime import datetime

def get_logger(name, log_dir='logs'):
    """
    콘솔과 파일에 동시에 로그를 남기는 로거를 생성하여 반환합니다.
    
    Args:
        name (str): 로거 이름 (보통 __name__ 사용)
        log_dir (str): 로그 파일을 저장할 디렉토리 경로
    
    Returns:
        logging.Logger: 설정된 로거 객체
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
    
    # 2. 파일 핸들러 (FileHandler)
    try:
        # 로그 디렉토리 생성
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
            
        # 날짜별 로그 파일명 생성 (예: logs/2024-05-20_analysis.log)
        today = datetime.now().strftime('%Y-%m-%d')
        log_filepath = os.path.join(log_dir, f"{today}_analysis.log")
        
        file_handler = logging.FileHandler(log_filepath, encoding='utf-8')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.INFO) # 파일에는 DEBUG 레벨까지 저장하려면 여기를 수정
        logger.addHandler(file_handler)
        
    except Exception as e:
        print(f"[Warning] 로그 파일 생성 실패: {e}")
    
    # 전파 방지 (Jupyter Notebook 중복 출력 방지)
    logger.propagate = False
    
    return logger