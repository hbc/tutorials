## Single-cell RNA-seq analysis tutorials

This repository contains tutorials for how to perform each part of a single-cell RNA-seq analysis, from running bcbio on the raw data to performing clustering, marker identificaton and differential expression analysis with DESeq2 and EdgeR. It also contains documents for generating data to use with online exploratory tools such as SPRING.

### Generating abundance estimates with bcbio

- [Running bcbio](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md)

### Setting up to run scRNA-seq analyses using R on O2

- [R on O2](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/R_set-up.md)

### Analysis workflow with Seurat (version 3)
repo in development at https://github.com/hbctraining/scRNA-seq/tree/master/lessons - will update links below when completed in mid-July
- [Quality control analysis]()
- [Clustering analysis]() 
- [Marker identification analysis]()
- [Integration analysis through marker identification]()

### Analysis workflow with Seurat (version 2)

- [Quality control analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_quality_control_analysis.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_QC_template.Rmd))
- [Clustering analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_clustering_analysis.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_clustering_template.Rmd))
- [Marker identification analysis](https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionIV/lessons/SC_marker_identification.html) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_marker_identification_template.Rmd))


### Downstream analyses
- Differential expression analysis with DESeq2
  - Preparation of data for DE analysis with DESeq2 (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_prep_for_DESeq2_analysis.Rmd))
  - DE analysis report (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/sc_DESeq2_analysis_report_template.Rmd))
- Differential expression analysis with EdgeR
- [Generating data for SPRING](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/SPRING.md)
- [Monocle trajectory analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/Monocle.md) (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/))
- Power analysis (associated `.Rmd` template available [here](https://github.com/hbc/tutorials/blob/master/scRNAseq/templates/))

### Analysis workflow with bcbioSingleCell [last update: 2017]

- [bcbioSingleCell set-up](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/bcbioSingleCell_setup.md)
- [Quality control analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/02_QC_report.md)
- [Clustering analysis](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/clustering_report_bcbioSingleCell.md)
