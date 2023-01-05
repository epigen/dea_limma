# Differential Analysis & Visualization Snakemake Workflow Using LIMMA
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for performing and visualizing differential expression (or accessibility) analyses (DEA) of NGS data (eg RNA-seq, ATAC-seq, scRNA-seq,...) powered by the R package [limma](https://www.bioconductor.org/packages/release/bioc/html/limma.html).

**If you use this workflow in a publication, don't forget to give credits to the authors by citing the URL of this (original) repository (and its DOI, see Zenodo badge above -> coming soon).**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

Table of contents
----------------
  * [Authors](#authors)
  * [Software](#software)
  * [Methods](#methods)
  * [Features](#features)
  * [Usage](#usage)
  * [Configuration](#configuration)
  * [Examples](#examples)
  * [Links](#links)
  * [Resources](#resources)

# Authors
- [Stephan Reichl](https://github.com/sreichl)
- [Lukas Folkman](https://github.com/lukas-folkman)


# Software
This project wouldn't be possible without the following software and their dependencies:

| Software       | Reference (DOI)                                   |
| :------------: | :-----------------------------------------------: |
| edgeR          | https://doi.org/10.1093/bioinformatics/btp616     |
| EnhancedVolcano| https://doi.org/10.18129/B9.bioc.EnhancedVolcano  |
| ggplot2        | https://ggplot2.tidyverse.org/                    |
| patchwork      | https://CRAN.R-project.org/package=patchwork      |
| pheatmap       | https://cran.r-project.org/package=pheatmap       |
| limma          | https://doi.org/10.1093/nar/gkv007                |
| Snakemake      | https://doi.org/10.12688/f1000research.29032.2    |

# Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References [ref] to the respective publications are curated in the software table above. Versions (ver) have to be read out from the respective conda environment specifications (.yaml file) or post execution. Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g. [X].

__Differential Expression Analysis (DEA).__ DEA was performed on the quality-controlled filtered [raw/normalized] counts using the LIMMA (ver) [ref] workflow for fitting a linear model [formula] to identify features (genes/regions) that statistically significantly change with [comparisons] compared to the control group [reference levels] (intercept). Briefly, we determined normalization factors with edgeR::calcNormFactors (optional) using method [X], then applied voom (optional) to estimate the mean-variance relationship of the log-counts. We used blocking on (optional) variable [X] to account for repeated measurements, lmFit to fit the model to the data, and finally eBayes (optional) with the robust (and trend flag – optional for normalized data) flag to compute (moderated/ordinary) t-statistics. For each comparison we extracted feature-wise average expression, effect sizes (log2 fold change) and their statistical significance as adjusted p-values, determined using the Benjamini-Hochberg method. Next, these results were filtered for relevant features based on the following criteria: statistical significance (adjusted p-value < [X]), absolute log2 fold change (> [X]), and average gene expression (> [X]). Finally, we performed hierarchical clustering on the effect sizes (log2 fold changes) of the union of all relevant features and comparison groups.

__Visualization.__ The filtered result statistics, i.e., number of relevant features split by positive (up) and negative (down) effect sizes, were visualized with stacked bar plots using ggplot (ver) [ref].
To visually summarize results of all performed comparisons, the filtered effect size (log2 fold change) values of all features that were found to be relevant in at least one comparison were plotted in a hierarchically clustered heatmap using pheatmap (ver) [ref]. 
Volcano plots were generated for each comparison using EnhancedVolcano (ver) [ref] with adjusted p-value threshold of [pCutoff] and log2 fold change threshold of [FCcutoff] as visual cut-offs for the y- and x-axis, respectively.
Finally, quality control plots of the fitted mean-variance relationship and raw p-values of the features were generated.


**The analysis and visualizations described here were performed using a publicly available Snakemake [ver] (ref) workflow [ref - cite this workflow here].**

# Features
The workflow performs the following steps that produce the outlined results:

- Differential Expression Analysis (DEA) steps:
  - (optional) calculation of normalization factors using __edgeR::calcNormFactors__.
  - (optional) calculation of precision weights to model the mean-variance relationship in order to make linear models "applicable" to count data (weighted least squares) using __voom__.
  - (optional) __block__ on a "group" factor in case you have repeated measurements (generalized least squares).
      - example use-case: you have N donors and T timepoints for each donor and want to model donor specific information like age, but still want to account for the variable __donor__ ie ~&nbsp;*timepoint*&nbsp;+&nbsp;*age*&nbsp;+&nbsp;*donor* is overdetermined hence the formula becomes ~&nbsp;*timepoint*&nbsp;+&nbsp;*age* and you "block" on __donor__.
  - fit linear models (ordinary least squares) with the design derived from the configured formula (expects "normal" data) using __lmFit__.
      - the fitted model object is saved (lmfit_object.rds) for alternative downstream analyses or manual inspection e.g., contrasts.
  - (optional) estimate variance "better" using __eBayes__, with the robustness flag (robust=TRUE), by looking across all genes (i.e. shrunk towards a common value) and compute moderated t-statistics.
      - (optional) eBayes with __limma-trend__ (trend=TRUE)
  - extract all statistics for variables of interest (=configured comparisons) using __topTable__ (eg coefficients/effect size, statistical significance,...).
  - save a feature list per comparison group and direction of change (up/down) for downstream analyses (eg enrichment analysis) (TXT).
    - (optional) annotated feature list with suffix "_annot" (TXT).
- DEA result statistics: total number of statistically significant features and split by positive (up) and negative (down) change (CSV).
- DEA result filtering of features (eg genes) by 
  - statistical significance (<= adjusted p-value: adj_pval)
  - effect size (>= absolute log 2 fold change: lfc)
  - average expression (>= ave_expr) in the data
- Log Fold Change (LFC) matrix of filtered features by comparison groups (CSV).
  - (optional) annotated LFC matrix with suffix "_annot" (CSV)
- Visualizations
  - filtered DEA result statistics ie number of features and direction (stacked bar plots)
  - volanco plot per comparison with configured cut-offs for statistical significance (pCutoff) and effect size (FCcutoff)
  - clustered heatmap of the LFC matrix
  - quality control plots
      - (optional) voom mean-variance trend
      - (optional) intermediate mean-variance trend, in case of blocking and vooming
      - post-fitting mean-variance trend
      - raw p-value distributions (to check for p-value inflation in relation to average expression)

FYI
- Colons (":") in variable/group names (eg interactions) are replaced with double-underscores ("\_\_") downstream.
- As of now we do not feature more complex contrast scenarios than are supported via topTable, but the fitted linear model is saved for downstream analyses.


# Usage
Here are some tips for the usage of this workflow:
- Limma usage and best practices are not explained. For deatiled documentation, tutorials and insctructions see [Resources](#resources). To understand the implementation at hand see [limma.R](./workflow/scripts/limma.R).
- Perform your first run(s) with loose filtering options/cut-offs and use the same for visualization to see if further filtering is even necessary or useful.
- Test minimal complex configurations on a subset of your data.
- Start with simple models and try to understand the results.

# Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# Examples
--- COMING SOON ---

# Links
- [GitHub Repository](https://github.com/epigen/dea_limma/)
- [GitHub Page](https://epigen.github.io/dea_limma/)
- [Zenodo Repository (coming soon)]()
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/dea_limma)

# Resources
- [Bioconductor - limma](http://bioconductor.org/packages/release/bioc/html/limma.html) includes a 150 page userguides
- [R Manual on Model Formulae](https://stat.ethz.ch/R-manual/R-patched/library/stats/html/formula.html)
- [Bioconductor - RNAseq123 - Workflow](https://bioconductor.org/packages/release/workflows/html/RNAseq123.html)
- limma workflow tutorial RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR
    - [notebook](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html)
    - [paper](https://f1000research.com/articles/5-1408/v3)
- A guide to creating design matrices for gene expression experiments
    - [notebook](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)
    - [paper](https://f1000research.com/articles/9-1444/v1)
- voom: precision weights unlock linear model analysis tools for RNA-seq read counts
    - [paper](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29)
- scRNA-seq DEA benchmark paper featuring limma-voom and limma-trend as valid/best methods
    - [paper: Bias, robustness and scalability in single-cell differential expression analysis](https://www.nature.com/articles/nmeth.4612)
    - [code](https://github.com/csoneson/conquer_comparison)
- scRNA-seq DEA benchmark featuring limma-voom and limma-trend emphasizing pseudo-bulking
    - [paper: Confronting false discoveries in single-cell differential expression](https://www.nature.com/articles/s41467-021-25960-2)
- alternative/complementary DEA method: Linear Mixed Models (LMM)
    - [variancePartition](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1323-z)
    - [dream](https://academic.oup.com/bioinformatics/article/37/2/192/5878955)
