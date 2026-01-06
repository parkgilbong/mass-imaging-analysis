# Snakemake Pipeline for Mass Imaging Analysis

Complete Snakemake workflow to replace the Jupyter notebook pipeline.

## Quick Start

### 1. Install/Update Environment

```bash
# Create or update conda environment with Snakemake
conda env create -f environment.yml

# Or update existing environment
conda env update -f environment.yml --prune

# Activate environment
conda activate mass-imaging-analysis
```

### 2. Run the Complete Pipeline

```bash
# Run all steps (parse → aggregate → analyze)
snakemake --cores 1

# Dry run to see what will be executed
snakemake --dry-run --printshellcmds

# Run with more cores (if applicable)
snakemake --cores 4
```

### 3. Run Individual Steps

```bash
# Step 1: Parse data only
snakemake parse_data --cores 1

# Step 2: Aggregate data only
snakemake aggregate_data --cores 1

# Step 3: Analyze data only
snakemake analyze_data --cores 1
```

---

## Pipeline Structure

The Snakefile defines three main rules corresponding to the original pipeline steps:

### Rule: `parse_data` (Step 1)
- **Input**: `.imzML` files from `data_dir`, `config.yaml`, m/z bins file
- **Output**: `*_mean_intensities.csv` files for each sample
- **Command**: `python src/main.py --config config/config.yaml`

### Rule: `aggregate_data` (Step 2)
- **Input**: All `*_mean_intensities.csv` files
- **Output**: `aggregated_mean_intensities.csv`
- **Command**: `python src/aggregate.py --config config/config.yaml`

### Rule: `analyze_data` (Step 3)
- **Input**: `aggregated_mean_intensities.csv`
- **Output**: Statistical results and plots
- **Command**: `python src/analysis.py --config config/config.yaml`

---

## Features

### ✅ Automatic Dependency Management
Snakemake automatically determines which steps need to be run based on file timestamps and dependencies.

### ✅ Variable Serial Sections Support
The Snakefile fully supports the variable serial sections feature:
```yaml
group_info:
  - name: "saline"
    n_per_group: 3
    num_serial: [2, 1, 3]  # Different per individual
```

### ✅ Logging
All pipeline steps log their output to `output/<output_dir>/logs/`:
- `parse_data.log`
- `aggregate_data.log`
- `analyze_data.log`

### ✅ Reproducibility
The entire pipeline is defined in a single file, making it easy to reproduce results.

---

## Advanced Usage

### Custom Config File

The repository includes `config/config.yaml` as a template. To use your own configuration:

1. **Create your custom config** (e.g., `config/config_myproject.yaml`)
2. **Run with custom config**:

```bash
# Use a different config file
snakemake --configfile config/config_myproject.yaml --cores 1
```

> [!NOTE]
> Files matching `config/config_*.yaml` are automatically ignored by git, so you can create multiple project-specific configs without tracking them in version control. Only `config/config.yaml` (the template) is tracked.

### Force Re-run

```bash
# Force re-run all steps
snakemake --forceall --cores 1

# Force re-run specific rule
snakemake --forcerun parse_data --cores 1
```

### Clean Output

```bash
# Clean all output files
snakemake clean
```

### Visualize DAG

```bash
# Generate workflow diagram
snakemake --dag | dot -Tpng > dag.png

# Or use rulegraph for simplified view
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

---

## Comparison: Jupyter vs Snakemake

| Feature | Jupyter Notebook | Snakemake |
|---------|-----------------|-----------|
| **Execution** | Manual cell-by-cell | Automatic dependency-based |
| **Reproducibility** | Requires careful cell order | Guaranteed by DAG |
| **Partial Re-run** | Manual tracking | Automatic based on timestamps |
| **Parallelization** | Not supported | Built-in with `--cores` |
| **Logging** | Mixed with output | Separate log files |
| **Integration** | Interactive only | CLI, HPC, cloud-ready |
| **Debugging** | Interactive | Via log files |

---

## Troubleshooting

### Issue: "Nothing to be done"
**Cause**: All output files are up-to-date.
**Solution**: Use `--forceall` to re-run or delete output files.

### Issue: Missing input files
**Cause**: `.imzML` files not found in `data_dir`.
**Solution**: Check `config.yaml` paths and file naming convention.

### Issue: Config validation error
**Cause**: `num_serial` list length doesn't match `n_per_group`.
**Solution**: Fix config.yaml to ensure consistency.

---

## Migration from Jupyter Notebook

If you were using `main.ipynb`, the Snakemake pipeline provides the same functionality:

**Old workflow (Jupyter):**
```python
# Cell 1
import src.main as main
main.main()

# Cell 2
import src.aggregate as aggregate
aggregate.main()

# Cell 3
import src.analysis as analysis
analysis.main()
```

**New workflow (Snakemake):**
```bash
snakemake --cores 1
```

Both produce identical results, but Snakemake offers better:
- Dependency tracking
- Reproducibility
- Integration with HPC/cloud systems
- Automatic partial re-runs

---

## Next Steps

1. **Test the pipeline**: Run `snakemake --dry-run` to verify setup
2. **Run full pipeline**: Execute `snakemake --cores 1`
3. **Check logs**: Review log files in `output/<output_dir>/logs/`
4. **Customize**: Modify `Snakefile` for your specific needs

For more information, see the [Snakemake documentation](https://snakemake.readthedocs.io/).
