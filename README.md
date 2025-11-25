# Mass Imaging Analysis Pipeline

MSI (Mass Spectrometry Imaging) 데이터를 처리하고 통계 분석을 수행하는 파이프라인입니다. `.imzML` 형식의 원본 데이터를 파싱하고, 그룹 간 통계 비교를 수행하며, 시각화 결과를 생성합니다.

## 목차
- [시스템 요구사항](#시스템-요구사항)
- [설치 방법](#설치-방법)
  - [1. Miniconda 설치](#1-miniconda-설치)
  - [2. 환경 설정](#2-환경-설정)
- [사용 방법](#사용-방법)
  - [Jupyter Notebook 실행](#jupyter-notebook-실행)
  - [main.ipynb 실행](#mainipynb-실행)
- [파이프라인 구조](#파이프라인-구조)
  - [Step 1: 데이터 파싱 (main.py)](#step-1-데이터-파싱-mainpy)
  - [Step 2: 데이터 집계 (aggregate.py)](#step-2-데이터-집계-aggregatepy)
  - [Step 3: 통계 분석 및 시각화 (analysis.py)](#step-3-통계-분석-및-시각화-analysispy)
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

Miniconda는 Python과 conda 패키지 매니저가 포함된 경량 배포판입니다.

#### Windows
1. [Miniconda 다운로드 페이지](https://docs.conda.io/en/latest/miniconda.html)에서 Windows용 설치 프로그램을 다운로드합니다.
2. 다운로드한 `.exe` 파일을 실행하여 설치를 진행합니다.
3. 설치 중 "Add Anaconda to my PATH environment variable" 옵션은 **체크하지 않는 것을 권장**합니다.
4. 설치 완료 후, **Anaconda Prompt**를 실행합니다.

#### macOS
1. [Miniconda 다운로드 페이지](https://docs.conda.io/en/latest/miniconda.html)에서 macOS용 설치 프로그램을 다운로드합니다.
   - Intel 프로세서: `Miniconda3-latest-MacOSX-x86_64.sh`
   - Apple Silicon (M1/M2/M3): `Miniconda3-latest-MacOSX-arm64.sh`
2. 터미널을 열고 다운로드한 파일을 실행합니다:
   ```bash
   bash Miniconda3-latest-MacOSX-*.sh
   ```
3. 설치 과정에서 라이선스 동의 후 설치 경로를 확인합니다 (기본값 사용 권장).
4. 설치 완료 후, 터미널을 재시작합니다.

#### Linux
1. 터미널에서 다음 명령어를 실행합니다:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
2. 설치 과정을 진행하고, 터미널을 재시작합니다.

**설치 확인:**
```bash
conda --version
```

### 2. 환경 설정

프로젝트에 필요한 Python 패키지들을 포함하는 conda 환경을 생성합니다.

1. **프로젝트 디렉토리로 이동:**
   ```bash
   cd /path/to/mass-imaging-analysis
   ```

2. **conda 환경 생성:**
   ```bash
   conda env create -f environment.yml
   ```
   
   이 명령어는 `environment.yml` 파일에 정의된 대로 `mass-imaging-analysis`라는 이름의 환경을 생성하고, 다음 패키지들을 설치합니다:
   - Python 3.11
   - pyimzml (imzML 파일 파싱)
   - pandas, numpy (데이터 처리)
   - pyyaml (설정 파일 파싱)
   - scipy, statsmodels (통계 분석)
   - seaborn, matplotlib (시각화)
   - jupyterlab (노트북 실행)

3. **환경 활성화:**
   ```bash
   conda activate mass-imaging-analysis
   ```

4. **설치 확인:**
   ```bash
   python --version
   conda list
   ```

---

## 사용 방법

### 데이터 준비

1. **데이터 디렉토리 구조:**
   ```
   mass-imaging-analysis/
    ├── config/
    │   ├── config.yaml
    │   └── Mass ranges of molecules.csv
    ├── data/
    │   └── 4_groups/
    │       ├── saline 1-1 cortex-total ion count.imzML
    │       ├── saline 1-1 cortex-total ion count.ibd
    │       └── ...
    ├── logs/                  <-- (자동 생성) 실행 로그 저장 폴더
    ├── output/
    │   └── 4_groups/          <-- 결과 파일 저장 폴더
    ├── src/
    │   ├── utils/
    │   │   └── logging_utils.py
    │   ├── main.py
    │   ├── aggregate.py
    │   ├── analysis.py
    │   └── parsing.py
    ├── main.ipynb
    └── environment.yml
   ```

2. **설정 파일 수정:**
   - `config/config.yaml` 파일을 열어 실험 설정에 맞게 수정합니다.
   - 그룹명, 샘플 수, ROI 정보 등을 확인하고 조정합니다.

### Jupyter Notebook 실행

1. **환경 활성화 (아직 활성화하지 않은 경우):**
   ```bash
   conda activate mass-imaging-analysis
   ```

2. **JupyterLab 실행:**
   ```bash
   jupyter lab
   ```
   또는 Jupyter Notebook을 선호하는 경우:
   ```bash
   jupyter notebook
   ```

3. **브라우저 자동 실행:**
   - 명령어 실행 후 자동으로 브라우저가 열립니다.
   - 만약 열리지 않으면, 터미널에 표시된 URL (예: `http://localhost:8888/...`)을 복사하여 브라우저에 붙여넣습니다.

### main.ipynb 실행

1. **JupyterLab/Notebook에서 `main.ipynb` 파일을 엽니다.**

2. **셀 순차 실행:**
   - 각 셀을 순서대로 실행합니다 (`Shift + Enter` 또는 상단의 "Run" 버튼).
   - 또는 전체 노트북을 한 번에 실행할 수 있습니다: `Cell` → `Run All`

3. **실행 순서:**
   - **셀 1-2**: 필요한 모듈 import 및 경로 설정
   - **셀 3**: `main.py` 실행 (데이터 파싱)
   - **셀 4**: `aggregate.py` 실행 (데이터 집계)
   - **셀 5**: `analysis.py` 실행 (통계 분석 및 시각화)
   - **셀 6-7**: 결과 확인 (CSV 파일 및 플롯 이미지)

4. **실행 시간:**
   - 파일 수와 크기에 따라 다르지만, 전체 파이프라인 실행에 수 분에서 수십 분이 소요될 수 있습니다.

### 로그 확인 
실행 중 발생하는 모든 메시지와 오류는 콘솔뿐만 아니라 파일로도 기록됩니다.

* 위치: logs/ 디렉토리
* 파일명: YYYY-MM-DD_analysis.log (예: 2024-05-21_analysis.log)
* 내용: 실행 시간, 처리 중인 파일, 경고(Warning) 및 오류(Error) 상세 내역
---

## 파이프라인 구조

### Step 1: 데이터 파싱 (main.py)

**목적:** `.imzML` 파일을 읽어 m/z bin별 intensity 데이터를 추출하고 CSV로 저장합니다.

**주요 기능:**
1. `config.yaml` 파일에서 실험 설정 로드
2. 처리할 모든 `.imzML` 파일 목록 생성 및 유효성 검사
3. m/z bin 정보 로드 (`Mass ranges of molecules.csv` 또는 config의 `direct_bins`)
4. 각 `.imzML` 파일에 대해:
   - 모든 스펙트럼을 읽어 m/z bin별로 intensity 합산
   - 전체 intensity 매트릭스를 CSV로 저장 (`*_binned_spectra.csv`)
   - m/z bin별 평균 intensity를 CSV로 저장 (`*_mean_intensities.csv`)

**입력 파일:**
- `data/{data_dir}/*.imzML` 및 `*.ibd` 파일
- `config/config.yaml`
- `config/Mass ranges of molecules.csv` (binning mode가 'file'인 경우)

**출력 파일:**
- `output/{output_dir}/*_binned_spectra.csv`: 각 픽셀의 m/z bin별 intensity 매트릭스
  - 컬럼: `x`, `y`, `bin1_name`, `bin2_name`, ...
- `output/{output_dir}/*_mean_intensities.csv`: m/z bin별 평균 intensity (1행)
  - 컬럼: `bin1_name`, `bin2_name`, ...

**실행 코드 (main.ipynb 내):**
```python
import src.main as main
main.main()
```

---

### Step 2: 데이터 집계 (aggregate.py)

**목적:** Technical replicates (num_serial)의 평균을 계산하여 biological replicates 단위로 데이터를 집계합니다.

**주요 기능:**
1. Step 1에서 생성된 `*_mean_intensities.csv` 파일들을 로드
2. 동일한 (group, n, roi) 조합에 대해 여러 serial (s=1, s=2, ...)의 평균 계산
3. 모든 biological replicates를 하나의 DataFrame으로 통합
4. 최종 집계 결과를 CSV로 저장

**입력 파일:**
- `output/{output_dir}/*_mean_intensities.csv` (Step 1에서 생성된 파일들)
- `config/config.yaml`

**출력 파일:**
- `output/{output_dir}/aggregated_mean_intensities.csv`: 집계된 데이터
  - 각 행: 하나의 biological replicate (group, n, roi)
  - 컬럼: `group`, `n`, `roi`, `bin1_name`, `bin2_name`, ...
  - 예시: 4개 그룹 × 3개 n × 1개 roi = 12행

**실행 코드 (main.ipynb 내):**
```python
import src.aggregate as aggregate
aggregate.main()
```

---

### Step 3: 통계 분석 및 시각화 (analysis.py)

**목적:** 그룹 간 통계 비교를 수행하고 결과를 시각화합니다.

**주요 기능:**
1. Step 2에서 생성된 `aggregated_mean_intensities.csv` 파일 로드
2. 데이터를 Long format으로 변환 (각 행이 하나의 측정값)
3. ROI 및 m/z bin별로 통계 검정 수행:
   - **2개 그룹:** t-test (parametric) 또는 Mann-Whitney U test (non-parametric)
   - **3개 이상 그룹:** ANOVA (parametric) 또는 Kruskal-Wallis test (non-parametric)
   - **Post-hoc:** Tukey HSD (parametric) 또는 Bonferroni-corrected Mann-Whitney U (non-parametric)
4. 각 ROI에 대해:
   - 막대 그래프 + 개별 데이터 포인트 시각화
   - 통계적으로 유의한 차이가 있는 경우 '*' 표시
   - PNG 이미지로 저장
5. GraphPad Prism 형식의 CSV 파일 생성

**입력 파일:**
- `output/{output_dir}/aggregated_mean_intensities.csv` (Step 2에서 생성)
- `config/config.yaml`

**출력 파일:**
- `output/{output_dir}/statistical_results_main.csv`: 메인 통계 검정 결과
  - 컬럼: `roi`, `m_z_bin`, `test_name`, `p_value`, `significant`
- `output/{output_dir}/statistical_results_posthoc.csv`: Post-hoc 검정 결과
  - 컬럼: `roi`, `m_z_bin`, `test_name`, `group1`, `group2`, `p_value`, `p_adj`, `significant`
- `output/{output_dir}/plot_roi_{roi_name}.png`: ROI별 시각화 결과
- `output/{output_dir}/aggregated_data_roi_{roi_name}_prism.csv`: GraphPad Prism용 데이터

**실행 코드 (main.ipynb 내):**
```python
import src.analysis as analysis
analysis.main()
```

---

## 설정 파일

### config/config.yaml

프로젝트의 주요 설정을 정의하는 파일입니다.

```yaml
# 데이터 및 출력 경로
data_dir: "data/4_groups"           # .imzML 파일들이 있는 디렉토리
output_dir: "output/4_groups"       # 결과 파일이 저장될 디렉토리

# 그룹 정보
group_info:
  num_of_group: 4
  name_of_group:
    - "saline"
    - "glyoxylate"
    - "acetate"
    - "etoh"

n_per_group: 3                     # 각 그룹의 biological replicates 수
num_serial: 2                      # 각 biological replicate의 technical replicates 수

# ROI 정보
roi_info:
  num_of_roi: 1
  name_of_roi:
    - "cortex"

# Binning 설정
binning_settings:
  mode: "file"                     # 'file' 또는 'direct'
  file_path: "config/Mass ranges of molecules.csv"  # mode='file'일 때
  # direct_bins:                   # mode='direct'일 때 (여기서는 사용 안 함)
  #   - mz: 72.9926
  #     width: 0.003
  #     name: "Bin 1 (72)"

# 통계 설정
statistics_settings:
  test_type: "non_parametric"      # 'parametric' 또는 'non_parametric'
  p_value_threshold: 0.05          # 유의수준
```

**파일명 규칙:**
- 파일명 템플릿: `{group} {n}-{s} {roi}-total ion count.imzML`
- 예시: `saline 1-1 cortex-total ion count.imzML`

### config/Mass ranges of molecules.csv

SCiLS Lab에서 export한 m/z bin 정보 파일입니다.

**파일 형식:**
```
# 주석 줄들...
m/z;Interval Width (+/- Da);Color;Name;Intensity [Regions]
72.9926;0.003;#b2df8a;;224.01341247559
102.055;0.003;#cab2d6;;4542.3041992188
...
```

- `m/z`: bin의 중심 m/z 값
- `Interval Width (+/- Da)`: bin의 반폭 (±)
- `Name`: bin의 이름 (비어있으면 m/z 값 사용)

---

## 트러블슈팅

### 1. 환경 생성 실패

**증상:**
```
ResolvePackageNotFound: ...
```

**해결책:**
- 인터넷 연결 확인
- conda 업데이트: `conda update conda`
- 채널 우선순위 확인: `conda config --show channels`

### 2. 파일을 찾을 수 없음 오류

**증상:**
```
오류: 파일을 찾을 수 없습니다. data/.../xxx.imzML
```

**해결책:**
- `.imzML`과 `.ibd` 파일이 모두 존재하는지 확인
- `config.yaml`의 `data_dir` 경로가 올바른지 확인
- 파일명이 템플릿과 일치하는지 확인: `{group} {n}-{s} {roi}-total ion count.imzML`

### 3. Jupyter Notebook이 실행되지 않음

**증상:**
```
jupyter: command not found
```

**해결책:**
- 환경이 활성화되었는지 확인: `conda activate mass-imaging-analysis`
- JupyterLab 재설치: `conda install -c conda-forge jupyterlab`

### 4. 메모리 부족 오류

**증상:**
```
MemoryError: Unable to allocate array
```

**해결책:**
- 더 작은 데이터셋으로 테스트
- 시스템 메모리 확인 및 다른 프로그램 종료
- 한 번에 처리하는 파일 수 줄이기 (config 수정)

### 5. 통계 결과가 비어있음

**증상:**
Post-hoc 통계 결과 파일이 생성되지 않거나 비어있음

**해결책:**
- 메인 통계 검정에서 유의한 차이가 없는 경우 정상입니다.
- `config.yaml`의 `p_value_threshold` 값 조정 (예: 0.05 → 0.1)
- 데이터 품질 및 그룹 간 차이 확인

---

## 라이선스 및 인용

이 프로젝트를 사용하여 연구 결과를 발표하는 경우, 적절한 인용을 부탁드립니다.

---

## 문의

문제가 발생하거나 질문이 있는 경우, GitHub Issues를 통해 문의해주세요.
