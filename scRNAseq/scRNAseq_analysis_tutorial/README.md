## Single-cell RNA-seq analysis tutorials

This repository contains tutorials for how to perform each part of a single-cell RNA-seq analysis, from running bcbio on the raw data to performing clustering, marker identificaton and differential expression analysis with DESeq2 and EdgeR. It also contains documents for generating data to use with online exploratory tools such as SPRING.

### Generating abundance estimates with bcbio

- [Running bcbio](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md)

### Setting up to run scRNA-seq analyses using R on O2

- [R on O2](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/R_set-up.md)

### Analysis workflow with Seurat (version 3)
All steps from QC to integration, clustering, and marker identification can be found in the teaching team repo: https://hbctraining.github.io/scRNA-seq/schedule/. The hands-on lessons from the workshop can be found below:

- [Generation of count matrix](https://hbctraining.github.io/scRNA-seq/lessons/02_SC_generation_of_count_matrix.html)
- [Quality control set-up](https://hbctraining.github.io/scRNA-seq/lessons/03_SC_quality_control-setup.html)
- [Quality control](https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html)
- [Normalization and Integration](https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html)
- [Clustering](https://hbctraining.github.io/scRNA-seq/lessons/07_SC_clustering_cells_SCT.html)
- [Clustering QC](https://hbctraining.github.io/scRNA-seq/lessons/08_SC_clustering_quality_control.html)
- [Marker Identification](https://hbctraining.github.io/scRNA-seq/lessons/09_merged_SC_marker_identification.html)

### Downstream analyses
- [Generating data for SPRING](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/SPRING.md)
- [Differential expression analysis - pseudobulk method with DESeq2](https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html)
- [Velocity analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/velocity.md)
- Trajectory analysis: Slingshot
  - [Analysis/benchmarking Rmd report](https://github.com/hbc/hbc_scrnaseq_tseng_10x_brown_fat_mouse_hbc03764/blob/master/2019_09_tseng_multisample_analysis/analysis_reports/slingshot/tseng_slingshot_comprehensive_report.tar.gz)
  - [R script using Seurat clusters and UMAP dim reduction](https://github.com/hbc/hbc_scrnaseq_tseng_10x_brown_fat_mouse_hbc03764/blob/master/2019_09_tseng_multisample_analysis/analysis_reports/slingshot/VSM_only_slingshot_to_adipo24_UMAP.R)
- Power analysis 
  - [Current] Associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/)

## Past analysis workflows (deprecated)

### Power analysis prior to pseudobulk technique

- [Deprecated] Preparation of data for DE analysis with DESeq2 (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_prep_for_DESeq2_analysis.Rmd)) - Needs to be changed to pseudobulk analysis
- [Deprecated] DE analysis report (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_DESeq2_analysis_report_template.Rmd)) - Needs to be changed to pseudobulk analysis
  
### Analysis workflow with Seurat (version 2)

- [Quality control analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_quality_control_analysis.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_QC_template.Rmd))
- [Clustering analysis](https://hbctraining.github.io/scRNA-seq/lessons/05_SC_clustering_cells.html)
- [Marker identification analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_marker_identification.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_marker_identification_template.Rmd))

### Analysis workflow with bcbioSingleCell [last update: 2017]

- [bcbioSingleCell set-up](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/bcbioSingleCell_setup.md)
- [Quality control analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/02_QC_report.md)
- [Clustering analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/clustering_report_bcbioSingleCell.md)
