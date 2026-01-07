# Mass Imaging Analysis Pipeline

MSI (Mass Spectrometry Imaging) 데이터를 처리하고 통계 분석을 수행하는 파이프라인입니다. `.imzML` 형식의 원본 데이터를 파싱하고, 그룹 간 통계 비교를 수행하며, 시각화 결과 및 HTML 리포트를 생성합니다.

**최신 기능:**
- **Snakemake 워크플로우 지원:** 전체 분석 과정을 자동화하고 병렬 처리를 지원합니다.
- **유연한 Serial ID:** 비순차적인 serial section 번호를 지원합니다.
- **HTML 리포트:** 분석 결과를 요약, 그래프, 통계 테이블이 포함된 대화형 HTML 리포트로 제공합니다.

## 목차
- [시스템 요구사항](#시스템-요구사항)
- [설치 방법](#설치-방법)
- [사용 방법 (Snakemake)](#사용-방법-snakemake-권장)
- [사용 방법 (Jupyter Notebook)](#사용-방법-jupyter-notebook)
- [주요 기능](#주요-기능)
- [파이프라인 구조](#파이프라인-구조)
- [설정 파일](#설정-파일)
- [트러블슈팅](#트러블슈팅)

---

## 시스템 요구사항

- **운영체제**: Windows, macOS, Linux
- **Python**: 3.11
- **메모리**: 최소 8GB RAM 권장
- **디스크 공간**: 데이터 크기에 따라 다름 (수 GB 이상)

---

## 설치 방법

### 1. Miniconda 설치
(Miniconda 또는 Anaconda가 이미 설치되어 있다면 건너뛰세요.)

### 2. 환경 설정

1. **프로젝트 디렉토리로 이동:**
   ```bash
   cd /path/to/mass-imaging-analysis
   ```

2. **conda 환경 생성:**
   ```bash
   conda env create -f environment.yml
   ```

3. **환경 활성화:**
   ```bash
   conda activate mass-imaging-analysis
   ```

---

## 사용 방법 (Snakemake, 권장)

Snakemake를 사용하면 전체 파이프라인을 효율적으로 실행하고 관리할 수 있습니다.

### 1. 설정 파일 수정
`config/config.yaml` 파일을 열어 실험 설정(데이터 경로, 그룹 정보, ROI 등)을 수정합니다.

### 2. 워크플로우 실행
터미널에서 다음 명령어를 실행합니다:

```bash
# Dry-run (실행 계획 확인)
snakemake --dry-run

# 전체 파이프라인 실행 (코어 4개 사용)
snakemake --cores 4
```

더 자세한 사용법은 [SNAKEMAKE_GUIDE.md](SNAKEMAKE_GUIDE.md)를 참고하세요.

---

## 사용 방법 (Jupyter Notebook)

대화형 분석을 원하시면 Jupyter Notebook을 사용할 수 있습니다.

1. **JupyterLab 실행:**
   ```bash
   jupyter lab
   ```

2. **`main.ipynb` 실행:**
   노트북의 셀을 순차적으로 실행하여 분석을 진행합니다.

---

## 주요 기능

### 1. 유연한 Serial ID 지원
개체별로 서로 다른 수의 serial section이나 비순차적인 ID(예: 1번 누락, 2번만 존재)를 처리할 수 있습니다.

```yaml
# config.yaml 예시
group_info:
  - name: "treatment"
    n_per_group: 3
    serial_ids:
      - [1, 2]    # Mouse 1
      - [2]       # Mouse 2 (1번 누락)
      - [1, 2, 3] # Mouse 3
```

### 2. HTML 분석 리포트
분석 완료 후 `output/{output_dir}/analysis_report.html` 파일이 생성됩니다.
- **Dashboard:** 분석 요약 정보
- **Plots:** ROI별 Montage Plot (확대 가능)
- **Tables:** 검색 및 정렬 가능한 통계 결과 테이블
- **Configuration:** 분석에 사용된 설정 정보

---

## 파이프라인 구조

### Step 1: 데이터 파싱 (`src/parse_data.py`)
`.imzML` 파일을 읽어 m/z bin별 intensity 데이터를 추출하고 CSV로 저장합니다.

### Step 2: 데이터 집계 (`src/aggregate_data.py`)
Technical replicates (serial sections)의 평균을 계산하여 Biological replicates 단위로 데이터를 집계합니다.

### Step 3: 통계 분석 (`src/analyze_stats.py`)
그룹 간 통계 비교(T-test, ANOVA 등)를 수행하고 결과를 CSV로 저장합니다.

### Step 4: 시각화 및 리포팅 (`src/visualize_data.py`)
통계 결과를 바탕으로 Bar plot을 생성하고, 최종 HTML 리포트를 작성합니다.

---

## 설정 파일 (`config/config.yaml`)

프로젝트의 모든 설정은 이 파일에서 관리합니다. 자세한 주석이 파일 내에 포함되어 있습니다.

---

## 트러블슈팅


---

## 문의
문제가 발생하거나 질문이 있는 경우, GitHub Issues를 통해 문의해주세요.
