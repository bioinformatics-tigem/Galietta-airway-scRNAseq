# Galietta airway scRNA-seq

This repository contains a Scanpy notebook for processing and exploring single-cell RNA-seq data from differentiated human bronchial epithelial cells (HBECs). The dataset includes non-cystic fibrosis (non-CF) and cystic fibrosis (CF) samples, as well as non-CF cells exposed to inflammatory stimuli.

The **full processed dataset**, including all four samples, is available from the [Gene Expression Omnibus (GEO), accession GSE335566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335566). The associated [TRPV4–CFTR study](https://doi.org/10.1016/j.isci.2026.117304) focuses on the untreated non-CF sample for its single-cell analysis, but the GEO deposit contains the other conditions too.

## Dataset overview

| Sample ID | Sample type | Condition | Cells after QC |
| --- | --- | --- | ---: |
| [`BE63_NON_CF_NT`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9815718) | Non-CF | Untreated | 14,086 |
| [`BE73_NT`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9815719) | CF | Untreated | 16,090 |
| [`BE63_NON_CF_TNFa_IL17`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9815720) | Non-CF | TNFα + IL-17 | 8,574 |
| [`BE63_NON_CF_IFNa`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9815721) | Non-CF | IFNα | 12,338 |
| **Total** | | | **51,088** |

The dataset contains **51,088 cells and 36,601 genes after quality control**. BE63 and BE73 are distinct biological samples; BE63 appears under three conditions, while BE73 is an untreated CF sample. These four entries are not four independent donor replicates.

The four sample-level 10x Genomics filtered count matrices are also available through their individual GEO sample records. Together they contain 58,883 cells before the additional QC filtering performed in Scanpy.

## What is in the processed AnnData file?

Download `GSE335566_adata_GEO.h5ad` from GEO. It contains the count data, sample information, cell-type annotations and the results of the main dimensionality-reduction and clustering steps.

| AnnData field | Contents |
| --- | --- |
| `.X` | Gene expression normalized to 10,000 counts per cell and transformed with `log1p` |
| `.layers["counts"]` | Original gene counts from the 10x matrices, retained after QC |
| `.raw.X` | Log-normalized expression values, **not** raw counts |
| `.obs["sample_id"]` | Sample identity for each cell |
| `.obs["curated_coarse_cell_type"]` | Broad epithelial cell-type annotations |
| `.obs["curated_cell_type"]` | More detailed cell-type annotations |
| `.obs["curated_cl_term"]`, `.obs["curated_cl_term_id"]` | Cell Ontology terms and identifiers |
| `.obs["leiden_nn30_res_0_5"]` | Original 13 Leiden clusters |
| `.var["highly_variable"]` | The 3,000 genes selected for dimensionality reduction |
| `.obsm["X_pca"]`, `.obsm["X_umap"]` | Precomputed PCA and UMAP coordinates |
| `.obsp["distances"]`, `.obsp["connectivities"]` | Precomputed nearest-neighbor graph |

All 36,601 genes are retained in `.X`. The 3,000 highly variable genes were used for PCA and clustering, but do not limit which genes can be inspected in the processed object.

**Important:** Use `.layers["counts"]` when you need the original gene counts. Despite its name, `.raw.X` stores log-normalized expression in this object.

The PCA, neighbor graph, UMAP and Leiden clusters were computed using the combined dataset. Selecting one sample keeps its coordinates and cluster labels from the combined analysis; it does not recompute the embedding.

## Getting started

Install Scanpy and download the processed H5AD from GEO. Then open it in Python:

```python
import scanpy as sc

adata = sc.read_h5ad("GSE335566_adata_GEO.h5ad")

print(adata)
print(adata.obs["sample_id"].value_counts())

# Explore the existing cell-type annotations
sc.pl.umap(
    adata,
    color="curated_coarse_cell_type",
    legend_loc="right margin",
    frameon=False,
)
```

To select the untreated non-CF sample used in the TRPV4–CFTR analysis:

```python
adata_nt = adata[adata.obs["sample_id"] == "BE63_NON_CF_NT"].copy()
```

## Analysis notebook

[`notebooks/SingleCell.ipynb`](notebooks/SingleCell.ipynb) documents the main workflow: reading the four 10x H5 files, QC, normalization, highly variable gene selection, PCA, nearest-neighbor graph construction, UMAP, Leiden clustering, manual cell-type annotation and example expression plots. The code and sample metadata are in `notebooks/` and `configs/`; package requirements are listed in [`requirements.txt`](requirements.txt).

If you want to **explore the published annotations and figures**, start from the processed H5AD: you do not need to rerun clustering. 

If you want to **rerun the pipeline from the four 10x H5 files**, place them in the configured data directory (`data/processedH5Files`) and check the filename pattern in `configs/config.yaml`. Files downloaded directly from GEO have a `GSM..._` prefix, so the input pattern must be adapted (or the files renamed) before running the notebook. Run Jupyter from the repository root so relative data paths resolve correctly. Recomputed graph-based clusters may vary slightly across software environments; the GEO H5AD preserves the original labels and coordinates.

## Publication and data access

- **Processed data and individual sample count matrices:** [GEO GSE335566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335566)
- **Associated publication:** [*TRPV4-driven ATP release activates CFTR through ADORA2B and P2RY2 receptors in the airway epithelium*](https://doi.org/10.1016/j.isci.2026.117304)
