# Visualization Notebooks

This directory contains Jupyter notebooks for developing and refining visualizations.

## Notebooks

### 01_volcano_plot_development.ipynb
Interactive development of volcano plots for differential metabolite analysis.

**Features:**
- Load integrated results
- Interactive parameter tuning (thresholds, colors, sizes)
- Multiple comparison support
- Export publication-ready figures

**Usage:**
```bash
jupyter notebook notebooks/01_volcano_plot_development.ipynb
```

### 02_heatmap_development.ipynb (Coming soon)
Interactive development of heatmaps for metabolite clustering.

---

## Workflow

1. **Develop**: Experiment with parameters in notebooks
2. **Refine**: Iterate until satisfied with visualization
3. **Modularize**: Extract to `src/visualize.py` when ready
4. **Integrate**: Add to pipeline if needed

---

## Requirements

Ensure you have Jupyter installed:
```bash
conda install jupyter
# or
pip install jupyter
```

Additional visualization packages:
```bash
conda install matplotlib seaborn
```
