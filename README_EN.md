# Mass Imaging Analysis Pipeline

A pipeline for processing and statistically analyzing MSI (Mass Spectrometry Imaging) data. It parses `.imzML` format raw data, performs statistical comparisons between groups, and generates visualization results.

## Table of Contents
- [System Requirements](#system-requirements)
- [Installation](#installation)
  - [1. Installing Miniconda](#1-installing-miniconda)
  - [2. Environment Setup](#2-environment-setup)
- [Usage](#usage)
  - [Running Jupyter Notebook](#running-jupyter-notebook)
  - [Command Line Interface (CLI)](#command-line-interface-cli)
- [Pipeline Structure](#pipeline-structure)
  - [Step 1: Data Parsing (main.py)](#step-1-data-parsing-mainpy)
  - [Step 2: Data Aggregation (aggregate.py)](#step-2-data-aggregation-aggregatepy)
  - [Step 3: Statistical Analysis and Visualization (analysis.py)](#step-3-statistical-analysis-and-visualization-analysispy)
- [Configuration Files](#configuration-files)
- [Troubleshooting](#troubleshooting)

---

## System Requirements

- **Operating System**: Windows, macOS, Linux
- **Python**: 3.11
- **Memory**: Minimum 8GB RAM recommended
- **Disk Space**: Depends on data size (several GB or more)

---

## Installation

### 1. Installing Miniconda

Miniconda is a lightweight distribution that includes Python and the conda package manager.

#### Windows
1. Download the Windows installer from the [Miniconda download page](https://docs.conda.io/en/latest/miniconda.html).
2. Run the downloaded `.exe` file to proceed with installation.
3. During installation, it is **recommended NOT to check** the "Add Anaconda to my PATH environment variable" option.
4. After installation, launch **Anaconda Prompt**.

#### macOS
1. Download the macOS installer from the [Miniconda download page](https://docs.conda.io/en/latest/miniconda.html).
   - Intel processors: `Miniconda3-latest-MacOSX-x86_64.sh`
   - Apple Silicon (M1/M2/M3): `Miniconda3-latest-MacOSX-arm64.sh`
2. Open Terminal and run the downloaded file:
   ```bash
   bash Miniconda3-latest-MacOSX-*.sh
   ```
3. Accept the license agreement and confirm the installation path (default recommended).
4. Restart Terminal after installation.

#### Linux
1. Run the following commands in Terminal:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
2. Complete the installation process and restart Terminal.

**Verify installation:**
```bash
conda --version
```

### 2. Environment Setup

Create a conda environment containing the required Python packages for the project.

1. **Navigate to project directory:**
   ```bash
   cd /path/to/mass-imaging-analysis
   ```

2. **Create conda environment:**
   ```bash
   conda env create -f environment.yml
   ```
   
   This command creates an environment named `mass-imaging-analysis` as defined in `environment.yml` and installs the following packages:
   - Python 3.11
   - pyimzml (imzML file parsing)
   - pandas, numpy (data processing)
   - pyyaml (configuration file parsing)
   - scipy, statsmodels (statistical analysis)
   - scikit-learn (dimensionality reduction: PCA, t-SNE)
   - seaborn, matplotlib (visualization)
   - jupyterlab, ipywidgets (notebook execution and interactive widgets)
   - umap-learn (UMAP dimensionality reduction)

3. **Activate environment:**
   ```bash
   conda activate mass-imaging-analysis
   ```

4. **Jupyter Lab extension setup (optional):**
   
   To enable interactive progress bars in Jupyter Lab, run:
   ```bash
   # Enable ipywidgets extension (automatic in Jupyter Lab 3.0+)
   jupyter labextension install @jupyter-widgets/jupyterlab-manager
   ```
   
   > **Note**: In Jupyter Lab 3.0 and above, ipywidgets are automatically enabled, so you can skip this step.

5. **Verify installation:**
   ```bash
   python --version
   conda list
   ```

---

## Usage

### Data Preparation

1. **Data directory structure:**
   ```
   mass-imaging-analysis/
   ├── config/
   │   ├── config.yaml
   │   └── Mass ranges of molecules.csv
   ├── data/
   │   └── 4_groups/
   │       ├── saline 1-1 cortex-total ion count.imzML
   │       ├── saline 1-1 cortex-total ion count.ibd
   │       ├── saline 1-2 cortex-total ion count.imzML
   │       ├── saline 1-2 cortex-total ion count.ibd
   │       └── ... (other .imzML and .ibd files)
   ├── output/
   │   └── 4_groups/
   ├── src/
   ├── main.ipynb
   └── environment.yml
   ```

2. **Modify configuration file:**
   - Open `config/config.yaml` and modify according to your experimental setup.
   - Verify and adjust group names, sample counts, ROI information, etc.

### Running Jupyter Notebook

1. **Activate environment (if not already activated):**
   ```bash
   conda activate mass-imaging-analysis
   ```

2. **Launch JupyterLab:**
   ```bash
   jupyter lab
   ```
   Or if you prefer Jupyter Notebook:
   ```bash
   jupyter notebook
   ```

3. **Browser auto-launch:**
   - The browser should open automatically after running the command.
   - If not, copy the URL displayed in the terminal (e.g., `http://localhost:8888/...`) and paste it into your browser.

### Running main.ipynb

1. **Open `main.ipynb` file in JupyterLab/Notebook.**

2. **Execute cells sequentially:**
   - Run each cell in order (`Shift + Enter` or "Run" button at the top).
   - Or run the entire notebook at once: `Cell` → `Run All`

3. **Execution order:**
   - **Cells 1-2**: Import required modules and set paths
   - **Cell 3**: Run `main.py` (data parsing)
   - **Cell 4**: Run `aggregate.py` (data aggregation)
   - **Cell 5**: Run `analysis.py` (statistical analysis and visualization)
   - **Cells 6-7**: Check results (CSV files and plot images)

4. **Execution time:**
   - Depending on file count and size, the entire pipeline may take several minutes to tens of minutes.

### Command Line Interface (CLI)

You can run the pipeline directly from the terminal instead of using Jupyter Notebook.

1. **Run individual steps:**
   ```bash
   # Activate environment
   conda activate mass-imaging-analysis
   
   # Step 1: Data parsing
   python src/main.py --config config/config.yaml
   
   # Step 2: Data aggregation
   python src/aggregate.py --config config/config.yaml
   
   # Step 3: Statistical analysis and visualization
   python src/analysis.py --config config/config.yaml
   ```

2. **View help:**
   ```bash
   python src/main.py --help
   python src/aggregate.py --help
   python src/analysis.py --help
   ```

3. **Use with Snakemake:**
   
   The CLI interface is well-suited for use with workflow management tools like Snakemake.
   
   ```python
   # Example Snakefile
   rule parse_data:
       input:
           config="config/config.yaml"
       output:
           directory("output/4_groups")
       shell:
           "python src/main.py --config {input.config}"
   
   rule aggregate_data:
       input:
           config="config/config.yaml",
           parsed="output/4_groups"
       output:
           "output/4_groups/aggregated_mean_intensities.csv"
       shell:
           "python src/aggregate.py --config {input.config}"
   
   rule analyze_data:
       input:
           config="config/config.yaml",
           aggregated="output/4_groups/aggregated_mean_intensities.csv"
       output:
           "output/4_groups/statistical_results_main.csv"
       shell:
           "python src/analysis.py --config {input.config}"
   ```

---

## Pipeline Structure

### Step 1: Data Parsing (main.py)

**Purpose:** Read `.imzML` files, extract intensity data for each m/z bin, and save as CSV.

**Key Functions:**
1. Load experimental settings from `config.yaml`
2. Generate list of all `.imzML` files to process and validate
3. Load m/z bin information (`Mass ranges of molecules.csv` or config's `direct_bins`)
4. For each `.imzML` file:
   - Read all spectra and sum intensities for each m/z bin
   - Save full intensity matrix as CSV (`*_binned_spectra.csv`)
   - Save average intensity per m/z bin as CSV (`*_mean_intensities.csv`)

**Input Files:**
- `data/{data_dir}/*.imzML` and `*.ibd` files
- `config/config.yaml`
- `config/Mass ranges of molecules.csv` (if binning mode is 'file')

**Output Files:**
- `output/{output_dir}/*_binned_spectra.csv`: Intensity matrix per pixel for each m/z bin
  - Columns: `x`, `y`, `bin1_name`, `bin2_name`, ...
- `output/{output_dir}/*_mean_intensities.csv`: Average intensity per m/z bin (1 row)
  - Columns: `bin1_name`, `bin2_name`, ...

**Execution code (in main.ipynb):**
```python
import src.main as main
main.main()
```

---

### Step 2: Data Aggregation (aggregate.py)

**Purpose:** Calculate averages of technical replicates (num_serial) and aggregate data by biological replicates.

**Key Functions:**
1. Load `*_mean_intensities.csv` files generated in Step 1
2. Calculate average across multiple serials (s=1, s=2, ...) for the same (group, n, roi) combination
3. Consolidate all biological replicates into a single DataFrame
4. Save final aggregated results as CSV

**Input Files:**
- `output/{output_dir}/*_mean_intensities.csv` (files generated in Step 1)
- `config/config.yaml`

**Output Files:**
- `output/{output_dir}/aggregated_mean_intensities.csv`: Aggregated data
  - Each row: one biological replicate (group, n, roi)
  - Columns: `group`, `n`, `roi`, `bin1_name`, `bin2_name`, ...
  - Example: 4 groups × 3 n × 1 roi = 12 rows

**Execution code (in main.ipynb):**
```python
import src.aggregate as aggregate
aggregate.main()
```

---

### Step 3: Statistical Analysis and Visualization (analysis.py)

**Purpose:** Perform statistical comparisons between groups and visualize results.

**Key Functions:**
1. Load `aggregated_mean_intensities.csv` file generated in Step 2
2. Convert data to long format (each row is one measurement)
3. Perform statistical tests for each ROI and m/z bin:
   - **2 groups:** t-test (parametric) or Mann-Whitney U test (non-parametric)
   - **3+ groups:** ANOVA (parametric) or Kruskal-Wallis test (non-parametric)
   - **Post-hoc:** Tukey HSD (parametric) or Bonferroni-corrected Mann-Whitney U (non-parametric)
4. For each ROI:
   - Visualize with bar graphs + individual data points
   - Mark statistically significant differences with '*'
   - Save as PNG image
5. Generate CSV files in GraphPad Prism format

**Input Files:**
- `output/{output_dir}/aggregated_mean_intensities.csv` (generated in Step 2)
- `config/config.yaml`

**Output Files:**
- `output/{output_dir}/statistical_results_main.csv`: Main statistical test results
  - Columns: `roi`, `m_z_bin`, `test_name`, `p_value`, `significant`
- `output/{output_dir}/statistical_results_posthoc.csv`: Post-hoc test results
  - Columns: `roi`, `m_z_bin`, `test_name`, `group1`, `group2`, `p_value`, `p_adj`, `significant`
- `output/{output_dir}/plot_roi_{roi_name}.png`: Visualization results per ROI
- `output/{output_dir}/aggregated_data_roi_{roi_name}_prism.csv`: Data for GraphPad Prism

**Execution code (in main.ipynb):**
```python
import src.analysis as analysis
analysis.main()
```

---

## Configuration Files

### config/config.yaml

File defining the main project settings.

```yaml
# Data and output paths
data_dir: "data/4_groups"           # Directory containing .imzML files
output_dir: "output/4_groups"       # Directory where results will be saved

# Group information
group_info:
  - name: "saline"
    n_per_group: 3
    num_serial: 2
  - name: "glyoxylate"
    n_per_group: 3
    num_serial: 2
  - name: "acetate"
    n_per_group: 3
    num_serial: 2
  - name: "etoh"
    n_per_group: 3
    num_serial: 2

# ROI information
roi_info:
  - name: "cortex"

# Binning settings
binning_settings:
  mode: "file"                     # 'file' or 'direct'
  file_path: "config/Mass ranges of molecules.csv"  # when mode='file'

# Statistics settings
statistics_settings:
  test_type: "non_parametric"      # 'parametric' or 'non_parametric'
  p_value_threshold: 0.05          # Significance level
```

**File naming convention:**
- Filename template: `{group} {n}-{s} {roi}-total ion count.imzML`
- Example: `saline 1-1 cortex-total ion count.imzML`

#### Variable Serial Sections per Individual

**Traditional approach (all individuals have same number of serial sections):**
```yaml
group_info:
  - name: "saline"
    n_per_group: 3
    num_serial: 2  # Integer: all mice have 2 sections
```

**New approach (different number of serial sections per individual):**
```yaml
group_info:
  - name: "saline"
    n_per_group: 3
    num_serial: [2, 1, 3]  # List: mouse 1=2, mouse 2=1, mouse 3=3 sections
```

**Mixed usage:**
```yaml
group_info:
  # Traditional approach
  - name: "control"
    n_per_group: 3
    num_serial: 2
  
  # New approach
  - name: "treatment"
    n_per_group: 3
    num_serial: [2, 1, 3]
```

**Validation rules:**
- If `num_serial` is a list, its length must match `n_per_group`
- All values must be positive integers
- Clear error messages are displayed if validation fails

### config/Mass ranges of molecules.csv

M/z bin information file exported from SCiLS Lab.

**File format:**
```
# Comment lines...
m/z;Interval Width (+/- Da);Color;Name;Intensity [Regions]
72.9926;0.003;#b2df8a;;224.01341247559
102.055;0.003;#cab2d6;;4542.3041992188
...
```

- `m/z`: Center m/z value of the bin
- `Interval Width (+/- Da)`: Half-width of the bin (±)
- `Name`: Name of the bin (uses m/z value if empty)

---

## Troubleshooting

### 1. Environment Creation Failed

**Symptom:**
```
ResolvePackageNotFound: ...
```

**Solution:**
- Check internet connection
- Update conda: `conda update conda`
- Check channel priority: `conda config --show channels`

### 2. File Not Found Error

**Symptom:**
```
Error: File not found. data/.../xxx.imzML
```

**Solution:**
- Verify both `.imzML` and `.ibd` files exist
- Check that `data_dir` path in `config.yaml` is correct
- Verify filename matches template: `{group} {n}-{s} {roi}-total ion count.imzML`

### 3. Jupyter Notebook Won't Run

**Symptom:**
```
jupyter: command not found
```

**Solution:**
- Verify environment is activated: `conda activate mass-imaging-analysis`
- Reinstall JupyterLab: `conda install -c conda-forge jupyterlab`

### 4. Out of Memory Error

**Symptom:**
```
MemoryError: Unable to allocate array
```

**Solution:**
- Test with smaller dataset
- Check system memory and close other programs
- Reduce number of files processed at once (modify config)

### 5. Empty Statistical Results

**Symptom:**
Post-hoc statistical results file is not generated or is empty

**Solution:**
- This is normal if there are no significant differences in the main statistical test
- Adjust `p_value_threshold` in `config.yaml` (e.g., 0.05 → 0.1)
- Check data quality and differences between groups

---

## License and Citation

If you publish research results using this project, please provide appropriate citation.

---

## Contact

For issues or questions, please contact us through GitHub Issues.
