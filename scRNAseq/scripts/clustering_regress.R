# Usage: This Rscript takes two arguments: the path to the pre-regressed rds file and the name of the interestingGroup sample used in the analysis.
# Rscript name_of_script "path/to/file.rds" "interestingGroup_name" 
# To change the variables to regress, alter line 23.

# Load libraries and provide data directory

library(Seurat)
library(tidyverse)
data_dir <- "data"
load(file.path(data_dir, "cycle.rda"))
set.seed(1454944673L)

# Use command line argument to specify sample extracted
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)

# Read in desired data

pre_regressed_seurat <- readRDS(args[1])

# Regress out the uninteresting sources of variation in the data

vars_to_regress <- c("nUMI", "S.Score", "G2M.Score")

seurat <- ScaleData(pre_regressed_seurat, vars.to.regress = vars_to_regress)

# Re-run the PCA plots and color by cell cycle phase

seurat <- RunPCA(
  seurat,
  pc.genes = c(s_genes, g2m_genes),
  do.print = FALSE)

PCAPlot(seurat, group.by= "Phase")

# Perform the scoring for all genes

seurat <- seurat %>%
  RunPCA(do.print = FALSE) %>%
  ProjectPCA(do.print = FALSE)

# Create elbow plot

PCElbowPlot(seurat)

# Determine the estimate for significant PCs

pct = seurat@dr$pca@sdev / sum(seurat@dr$pca@sdev) * 100
cum = cumsum(pct)
co1 = which(cum > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),
           decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
pcs = min(co1, co2) # change to any other number

# Find cell clusters for different resolutions

seurat <- FindClusters(
  seurat,
  dims.use = 1:pcs,
  force.recalc = TRUE,
  print.output = TRUE,
  resolution = c(0.1, 0.6, 0.8, 1.0, 1.2, 1.8)
  save.SNN = TRUE)

PrintFindClustersParams(seurat)

# Save clustered cells

saveRDS(seurat, file = file.path(data_dir, paste0("seurat_", args[2], "_", args[3], "_tsne.rds")))
                                                                                                                                                                                                                             1,1           All

