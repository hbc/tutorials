# Usage: this Rscript is using Seurat to subset data for a particular sample given as a command line argument, then will perform normalization, calculation of variable genes and significant PCs, finally it will score cells for cell cycle. The output is ready for the regress.R script. 

# The script expects as input the name of the interestingGroup used to subset samples. If prefer to use a different column in metadata to subset, can change line 22.

# To run:  Rscript name_of_script "interestingGroup_name"

library(Seurat)
library(tidyverse)
data_dir <- "data"
load(file.path(data_dir, "cycle.rda"))
set.seed(1454944673L)

# Load pre_regressed Seurat object                                                                                                
pre_regressed_seurat <- readRDS(file.path(data_dir, "seurat_raw.rds"))

# Use command line argument to specify sample extracted
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# Subset by sample using commandline argument - uncomment code below if desired
# pre_regressed_seurat <- SubsetData(pre_regressed_seurat,
                                cells.use = rownames(pre_regressed_seurat@meta.data)[which(pre_regressed_seurat@meta.data$interestingGroups == args[1])])


# Normalize counts for total cell expression and take log value                            
pre_regressed_seurat <- pre_regressed_seurat %>% NormalizeData(
  normalization.method = "LogNormalize",
  scale.factor = 10000)

# Find variable genes based on the mean-dispersion relationship based on z-score for dispersion. Recommended to set parameters as to mark visual outliers on dispersion plot - default for ~2,000 variable genes

pre_regressed_seurat =  pre_regressed_seurat %>%
  FindVariableGenes(
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    do.plot = FALSE)


pre_regressed_seurat = pre_regressed_seurat %>%
  ScaleData(model.use = "linear")


# Check number of variable genes to determine if correct parameters used  

length(x = pre_regressed_seurat@var.genes)


pre_regressed_seurat <- CellCycleScoring(
  pre_regressed_seurat,
  g2m.genes = g2m_genes,
  s.genes = s_genes)

pre_regressed_seurat = RunPCA(
  pre_regressed_seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)

saveRDS(pre_regressed_seurat, file = file.path(data_dir, paste0("seurat_pre_regress_", args[1], ".rds")))
