# Transcriptomic Analysis of *Opaque-2* Loss-of-Function in *Zea mays*

**Author:** Ben Martin
**Data source:** [EMBL-EBI Expression Atlas — E-CURD-41](https://www.ebi.ac.uk/gxa/experiments/E-CURD-41)  
**License:** Unlicensed — all rights reserved. Contact the author for permission to reuse.

---

## Table of Contents

- [Project Overview](#project-overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Data](#data)
- [Methods](#methods)
- [Results](#results)
- [Limitations](#limitations)
- [References](#references)

---

## Project Overview

*Opaque-2* (*O2*) is a bZIP transcription factor in *Zea mays* (maize) that, prior to a landmark 2018 study, was well known but poorly characterised. While *O2* had long been associated with the regulation of seed storage protein gene expression, earlier transcriptome analyses of *o2* endosperm suggested a far broader regulatory role — including processes involved in the synthesis and metabolism of carbohydrates and lipids. However, it remained unclear whether *O2* was directly or indirectly responsible for these wider effects.

This project uses publicly available RNA-seq data to perform a differential expression analysis comparing wildtype *Zea mays* endosperm with *Opaque-2* loss-of-function mutants, with the aim of characterising the transcriptional consequences of *O2* knockout. Gene Ontology (GO) term enrichment analysis is used to assess whether the differentially expressed genes cluster into biologically meaningful functional categories, helping to distinguish direct regulatory targets from downstream or indirect effects.

---

## Repository Structure

```
.
├── README.md                  # This document
├── analysis.R                 # Full analysis pipeline (end-to-end)
├── sessionInfo.txt            # R session and package version information
├── data/
│   ├── E-CURD-41-experiment-design.tsv # Experimental design (metadata)
│    ├── E-CURD-41-raw-counts.tsv        # Raw data generated experimentally (RNA-Seq counts)
│    └── genes.txt                       # MaizeGDB gene annotation file      
├── plots/
│   ├── GO_barplot_frequency.png
│   ├── GO_dotplot_directional.png
│   ├── GO_dotplot_ORA.png
│   ├── MA_plot.png
│   ├── volcano_plot.png
│   ├── PCA_plot.png
│   └── heatmap_top20.png
└── results/
    ├── GO_term_frequency_table.csv
    └── GO_ORA_results.csv
```

---

## Dependencies

This analysis was conducted in R. The following packages are required:

| Package | Purpose |
|---|---|
| `DESeq2` | Differential expression analysis |
| `ggplot2` | General visualisation |
| `pheatmap` | Heatmap generation |
| `EnhancedVolcano` | Volcano plot |
| `stringr` | String manipulation for GO term parsing |
| `dplyr` | Data manipulation |
| `tidyr` | Data reshaping |


For exact package versions used in this analysis, see [`sessionInfo.txt`](sessionInfo.txt).

To install all required packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "EnhancedVolcano"))

install.packages(c("ggplot2", "pheatmap", "stringr", "dplyr", "tidyr"))
```

---

## Data

### RNA-seq Expression Data

Raw count data for 6 *Zea mays* endosperm samples (3 wildtype, 3 *Opaque-2* loss-of-function mutants) were retrieved from the EMBL-EBI Expression Atlas:

> **Accession:** [E-CURD-41](https://www.ebi.ac.uk/gxa/experiments/E-CURD-41)

Reads were aligned to the *Zea mays* reference genome using **TopHat v2.0.9**.

### Gene Annotation

Gene annotations were obtained from the [MaizeGDB](https://www.maizegdb.org/) project. The annotation file (`data/genes.txt`) contains gene identifiers in the `Zm00001ebxxxxxx` format, alternative gene names, full gene descriptions, and associated GO terms.

---

## Methods

### Differential Expression Analysis

A `DESeqDataSet` object was constructed from raw count data using the design formula:

```r
design = ~ genotype
```

Variance stabilising transformation was applied using `vst(dds, blind = FALSE)` for use in exploratory visualisation (PCA and heatmap). Differential expression was assessed using DESeq2's default Wald test. Genes were considered significantly differentially expressed at an adjusted p-value threshold of **padj < 0.05** (Benjamini-Hochberg correction).

### GO Term Enrichment Analysis

GO terms were parsed from the MaizeGDB annotation file. Each GO term entry follows the format `GO:XXXXXXX=description`; terms were extracted using this structure as an anchor to avoid misparsing caused by commas embedded within GO term descriptions. Duplicate gene–term pairs were removed prior to analysis.

**Over-Representation Analysis (ORA)** was performed using a one-sided Fisher's exact test for each GO term, comparing its frequency among significantly DE genes against its frequency across all annotated genes as background. Resulting p-values were corrected using the Benjamini-Hochberg method. GO terms were considered significantly enriched at **padj < 0.05**.

### Visualisation

The following plots were produced:

- **PCA plot** — principal component analysis of VST-normalised counts, used to assess replicate clustering and confirm genotype as the dominant source of variance
- **MA plot** — log-ratio vs mean expression, used to assess the overall distribution and magnitude of fold changes
- **Volcano plot** — log2 fold change vs -log10(p-value), highlighting significantly DE genes with gene name labels
- **Heatmap** — row-scaled VST expression values for the 20 most significantly DE genes (by adjusted p-value), with hierarchical clustering of both genes and samples
- **GO term bar chart** — frequency of GO terms among DE genes, split by direction of regulation
- **GO ORA dot plot** — enriched GO terms plotted by odds ratio and significance
- **Directional GO dot plot** — enriched GO terms separated by up- and downregulation

---

## Results

### Differential Expression

- Total significantly DE genes (padj < 0.05): 1603
- Upregulated in *o2* mutant: 632
- Downregulated in *o2* mutant: 971

### Principal Component Analysis

The PCA plot of VST-normalised counts show that:
- Wildtype and *o2* mutant samples separate along PC1
- Replicates do not cluster well along PC2, with both wildtype and *o2* mutant samples showing one outlier along this axis

### GO Term Enrichment

- Significantly enriched GO terms (ORA padj < 0.05): 8
- Notable enriched biological processes (by Odds Ratio) included: nutrient reservoir activity, metal ion binding, valine biosynthetic process, carbohydrate synthetic process, acetolactate synthase activity, cytoplasm, protein storage vacuole, and plasma membrane

Full results are available in [`results/GO_ORA_results.csv`](results/GO_ORA_results.csv).

## GO Term Analysis
- Several GO terms (434) showed large effect sizes (odds ratio > 10) while not meeting the adjusted p-value threshold (p < 0.05). Most of these have very few gene counts, only 26 exist in this category with more than 2 gene counts, but these are terms that may be limited by statistical power arising from a small sample size. These may warrant investigation in future studies with greater replication.
- Of the GO terms that are significantly enriched, nutrient reservoir activity, carbohydrate metabolic processing, protein storage vacuole are roles that would be expected given its characterised role in regulation of zein genes. This is visible in 4 of the top 10 (sorted by adjusted p-value) most differential expressed genes.
- The other GO terms that are significantly enriched with large ORs include acetolactate synthase activity and valine biosynthetic process. These are both GOs with small numbers of genes ( < 5) which may be contributing to their larger ORs, but merit further investigation in a larger sample. 
- The terms with significant adjusted p-values and smaller ORs include cytoplasm, plasma membrane, and metal ion binding. The smaller ORs indicate that the GO terms are not entirely likely to appear in the DE genes compared to the background, and may be reflective of either the breadth of the annotation rather than anything specific about O2 regulation, or the process is moderately associated with O2 regulation rather than being a core or specific target, and their association should be interpreted with caution.
- While the results reported here differ from those reported by Zhan et al. (i.e. 1863 vs 1603 identified differentially expressed genes), they are broadly consistent. The differences could be caused by differences in DESeq2 version, low-count filtering thresholds, and any consequential effects on BH correction across the full gene set.
---

## Limitations

- **Low replicate number:** With 3 biological replicates per group, statistical power is limited. Fold change estimates for lowly expressed genes should be interpreted with caution.
- **Indirect effects:** This analysis identifies transcriptional differences between genotypes but cannot distinguish direct O2 targets from indirect downstream effects without complementary ChIP-seq or motif analysis data.
- **Annotation completeness:** Not all *Zea mays* genes in the MaizeGDB annotation carry GO term assignments. Genes lacking annotation are excluded from the ORA background, which may introduce a modest bias toward well-characterised gene families.
- **Single tissue/timepoint:** Data are derived from endosperm tissue at a single developmental stage. The regulatory role of O2 may differ in other tissues or developmental contexts.

--- 

## Extensions

- Currently, the GO term ORA includes only binary membership in terms of whether or not a gene is up- or down-regulated. Given the biological debate centers on whether or not O2 is a direct or indirect regulator of various cellular processes, the magnitude by which expression of genes within GO terms changes in wildtype vs mutant would be a worthwhile question to investigate. This would allow for distinction between a GO term where member genes are consistently strongly up- or down-regulated and one where change is weak or in inconsistent directionality.

--- 

## References

- Love MI, Huber W, Anders S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
- Trapnell C, et al. (2012). Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. *Nature Protocols*, 7, 562–578.
- EMBL-EBI Expression Atlas. Experiment E-CURD-41. https://www.ebi.ac.uk/gxa/experiments/E-CURD-41
- MaizeGDB. https://www.maizegdb.org/
- (Zhan J, Li G, Ryu CH, Ma C, Zhang S et al. (2018) Opaque-2 Regulates a Complex Gene Network Associated with Cell Differentiation and Storage Functions of Maize Endosperm.)
