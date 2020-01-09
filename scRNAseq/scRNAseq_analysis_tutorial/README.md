## Single-cell RNA-seq analysis tutorials

This repository contains tutorials for how to perform each part of a single-cell RNA-seq analysis, from running bcbio on the raw data to performing clustering, marker identificaton and differential expression analysis with DESeq2 and EdgeR. It also contains documents for generating data to use with online exploratory tools such as SPRING.

### Generating abundance estimates with bcbio

- [Running bcbio](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md)

### Setting up to run scRNA-seq analyses using R on O2

- [R on O2](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/R_set-up.md)

### Analysis workflow with Seurat (version 3)
All steps can be found in the teaching team repo: https://hbctraining.github.io/scRNA-seq/schedule/.
- [Quality control analysis](https://hbctraining.github.io/scRNA-seq/lessons/03_SC_quality_control.html)
- [Clustering analysis](https://hbctraining.github.io/scRNA-seq/lessons/05_SC_clustering_cells.html) 
- [Cluster exploration](https://hbctraining.github.io/scRNA-seq/lessons/06_SC_clustering_quality_control.html)
- [Marker identification analysis](https://hbctraining.github.io/scRNA-seq/lessons/07_SC_marker_identification.html)
- [Integration - clustering](https://hbctraining.github.io/scRNA-seq/lessons/08_SC_clustering_analysis_integration.html)
- [Integration - marker id](https://hbctraining.github.io/scRNA-seq/lessons/09_SC_marker_identification_integration.html)

### Analysis workflow with Seurat (version 2)

- [Quality control analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_quality_control_analysis.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_QC_template.Rmd))
- [Clustering analysis](https://hbctraining.github.io/scRNA-seq/lessons/05_SC_clustering_cells.html)
- [Marker identification analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_marker_identification.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_marker_identification_template.Rmd))


### Downstream analyses
- [Generating data for SPRING](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/SPRING.md)
- Differential expression analysis
  - Pseudobulk with EdgeR - see Meeta/Victor/Rory
  - Using MAST - see Victor
- Trajectory analysis
  - Slingshot - see Rory
- Power analysis (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/))
  - [Deprecated] Preparation of data for DE analysis with DESeq2 (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_prep_for_DESeq2_analysis.Rmd))
  - [Deprecated] DE analysis report (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_DESeq2_analysis_report_template.Rmd))

### Analysis workflow with bcbioSingleCell [last update: 2017]

- [bcbioSingleCell set-up](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/bcbioSingleCell_setup.md)
- [Quality control analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/02_QC_report.md)
- [Clustering analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/clustering_report_bcbioSingleCell.md)
