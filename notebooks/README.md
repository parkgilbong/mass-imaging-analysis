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

### 02_heatmap_development.ipynb
Interactive development of heatmaps for metabolite clustering and pattern visualization.

**Features:**
- Load integrated results data
- Multiple normalization methods (Z-score, min-max, robust scaling)
- Hierarchical clustering with customizable algorithms
- Dendrogram visualization (rows and columns)
- Interactive parameter tuning
- Data filtering (significance, fold change, variance)
- Export processed data for 3rd party tools
- Export publication-ready figures (PNG, PDF, SVG)

**Usage:**
```bash
jupyter notebook notebooks/02_heatmap_development.ipynb
```

### 03_dimensionality_reduction.ipynb
Interactive dimensionality reduction analysis for exploring metabolite patterns.

**Features:**
- Multiple algorithms: PCA, t-SNE, UMAP
- Normalization and log transformation options
- Interactive parameter tuning
- 2D scatter plot visualization with group coloring
- Variance explained plots (PCA)
- Algorithm comparison (side-by-side plots)
- Export plots and transformed data
- Comprehensive parameter guidelines

**Usage:**
```bash
jupyter notebook notebooks/03_dimensionality_reduction.ipynb
```

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
