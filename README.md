# Cellalyze

Cellalyze is a learning project for exploring single-cell RNA sequencing data in a browser. It uses Streamlit for the interface and Scanpy for basic single-cell analysis.

## Current Features

For AnnData `.h5ad` files, the app currently:

- Loads an uploaded dataset.
- Reports basic dataset information, including cell and gene counts.
- Normalizes counts to a target sum of 10,000 per cell.
- Applies `log1p` transformation.
- Selects up to 2,000 highly variable genes.
- Computes PCA and a nearest-neighbor graph.
- Attempts Leiden clustering and stores the result as `clusters`.
- Computes a UMAP embedding.
- Displays interactive Plotly UMAP visualizations.
- Colors UMAP cells by gene expression or cell metadata.
- Displays a Scanpy UMAP plot for the selected feature.

CSV and CSV.GZ files are currently supported as tabular previews. They can be inspected and plotted by column, but they are not yet converted into AnnData or processed through the single-cell workflow.

## Tech Stack

- Python
- Streamlit
- Scanpy and AnnData
- pandas
- NumPy
- Matplotlib
- Plotly

## Setup

Use Python 3.10 or newer. From the project directory, create a virtual environment and install the dependencies:

```bash
cd /Users/chloeyfang/scRNAseq-visualizer
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install streamlit scanpy pandas numpy matplotlib plotly
```

The existing `venv312` directory may point to a Python installation that is no longer present. Creating a fresh `venv` as shown above avoids that broken interpreter reference.

## Run the App

From the project directory, with the virtual environment activated:

```bash
streamlit run app.py
```

Streamlit will print a local URL, normally:

```text
http://localhost:8501
```

Stop the app with `Ctrl+C` in the terminal.

## Using the App

1. Open the local Streamlit URL in a browser.
2. Upload an `.h5ad` file from the sidebar.
3. Wait for preprocessing and UMAP calculation to finish.
4. Choose whether to color the UMAP by gene expression or cell metadata.
5. Select a gene or metadata column to explore.
6. Inspect the interactive plot and dataset summary.

A useful starter dataset should contain:

- An expression matrix in `adata.X`.
- Gene names in `adata.var_names`.
- Cell names in `adata.obs_names`.
- Optional cell annotations in `adata.obs`.

## Analysis Pipeline

The current H5AD workflow is:

```text
uploaded H5AD
    -> normalize total counts
    -> log1p transform
    -> highly variable gene selection
    -> PCA
    -> neighbor graph
    -> Leiden clusters
    -> UMAP
    -> interactive visualization
```

The original uploaded object is kept separate from the processed copy so the displayed raw-data summary is not replaced by the normalized values.

## Project Structure

```text
scRNAseq-visualizer/
├── app.py       # Streamlit application
└── README.md    # Project documentation
```

## Limitations and Next Steps

The project is intentionally early-stage. Good next improvements would be:

- Add marker-gene analysis with `sc.tl.rank_genes_groups`.
- Add quality-control metrics and plots, such as total counts, detected genes, and mitochondrial percentage.
- Add filtering controls before normalization.
- Convert count-matrix CSV files into AnnData with clear row and column orientation options.
- Add tabs for UMAP, clusters, marker genes, quality control, and dataset summary.
- Add controls for the number of highly variable genes, neighbors, principal components, and Leiden resolution.
- Add error handling for malformed files, datasets that are too small for PCA, and missing clustering dependencies.
- Add sample data and automated tests for the preprocessing functions.

## Learning Resources

- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [AnnData documentation](https://anndata.readthedocs.io/)
- [Streamlit documentation](https://docs.streamlit.io/)
# cellalyze
# cellalyze
