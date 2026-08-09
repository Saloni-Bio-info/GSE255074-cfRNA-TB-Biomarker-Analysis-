# Differential Gene Expression Analysis — GSE255074

Analysis of RNA-seq data from **GSE255074** to identify genes differentially expressed between tuberculosis (TB) positive and TB negative patients, followed by Gene Set Enrichment Analysis (GSEA). This project was completed during a summer internship at **IIT Hyderabad**.

## Overview

The pipeline takes raw RNA-seq count data and public GEO metadata, performs quality filtering and normalization, fits a linear model to identify differentially expressed genes (DEGs), and runs downstream functional enrichment (GO: Biological Process) on the ranked gene list.

## Dataset

- **Accession:** [GSE255074](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255074) (Gene Expression Omnibus)
- **Comparison:** TB positive vs. TB negative patients
- **Data type:** Bulk RNA-seq raw counts

## Workflow

1. **Data import**
   - Load raw count matrix (genes × samples) via `read.delim()`
   - Fetch sample metadata directly from GEO using `GEOquery::getGEO()`

2. **Preprocessing**
   - Strip Ensembl gene ID version suffixes
   - Remove low-count genes (`rowSums(counts) > 10`)
   - Align sample order between count matrix and metadata
   - Define condition groups (TB positive / TB negative) from metadata

3. **Filtering and normalization**
   - Build a `DGEList` object (`edgeR`)
   - Filter lowly expressed genes with `filterByExpr()`
   - Apply TMM normalization via `calcNormFactors()`
   - Visual QC with boxplots of logCPM before and after normalization

4. **Differential expression analysis**
   - Design matrix for the TB-positive vs. TB-negative contrast
   - `voomWithQualityWeights()` transformation to model mean-variance trend and sample-specific reliability
   - Linear modeling with `lmFit()`, contrast extraction with `contrasts.fit()`, and moderated statistics via `eBayes()` (`limma`)
   - DEGs called at `|log2FC| > 1` and `adjusted p-value < 0.05` (Benjamini-Hochberg)

5. **Output tables**
   - `All_DEGs.csv` — full results table
   - `Significant_DEGs.csv` — genes passing significance thresholds
   - `Upregulated_Genes.csv` / `Downregulated_Genes.csv` — split by direction of change

6. **Visualization**
   - Volcano plot of log2FC vs. −log10(adjusted p-value), colored by significance (`Volcano_Plot.png`)

7. **Gene Set Enrichment Analysis (GSEA)**
   - Convert Ensembl IDs to Entrez IDs (`org.Hs.eg.db`)
   - Rank genes by log2FC
   - Run `gseGO()` (GO Biological Process) via `clusterProfiler`
   - Dot plot of top enriched terms (`GSEA_DotPlot.png`)
   - Results table (`GSEA_GO_BP_Results.csv`)

## Requirements

**R packages:**
```r
install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "limma",
  "edgeR",
  "GEOquery",
  "biomaRt",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot"
))

install.packages("dplyr")
```

## Usage

1. Download the raw counts file for GSE255074 from GEO.
2. Open `DGE_Analysis1.R` in RStudio (or run via Rscript).
3. When prompted (`file.choose()`), select the downloaded counts file.
4. Run the script top to bottom — it will fetch metadata automatically, perform the analysis, and write all output files (CSVs and plots) to the working directory.

## Outputs

| File | Description |
|---|---|
| `All_DEGs.csv` | Complete differential expression results for all tested genes |
| `Significant_DEGs.csv` | DEGs meeting significance thresholds |
| `Upregulated_Genes.csv` | Significantly upregulated genes in TB-positive samples |
| `Downregulated_Genes.csv` | Significantly downregulated genes in TB-positive samples |
| `Volcano_Plot.png` | Volcano plot of all genes |
| `GSEA_GO_BP_Results.csv` | GO Biological Process enrichment results |
| `GSEA_DotPlot.png` | Dot plot of top enriched GO terms |

## Tools & Packages

- **limma** — linear modeling and differential expression testing
- **edgeR** — count filtering and TMM normalization
- **GEOquery** — retrieval of GEO sample metadata
- **biomaRt** — gene ID annotation
- **clusterProfiler** / **org.Hs.eg.db** / **enrichplot** — GSEA and visualization
- **ggplot2** / **dplyr** — plotting and data wrangling

## Acknowledgements

This analysis was carried out as part of a summer internship at the **Indian Institute of Technology (IIT) Hyderabad**.
