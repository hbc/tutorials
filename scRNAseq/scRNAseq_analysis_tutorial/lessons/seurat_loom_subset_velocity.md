# Subset loom file to cells of interest

1. Split Seurat file into separate groups if velocity analysis is split by group

```r
# Split Seurat object
obj.list <- SplitObject(seurat, 
                        split.by = "sample_simple")               
```

2. Save as separate objects

```r
# Extract separate objects and save associated cell IDs
ctrl <- obj.list[["ctrl"]]
ctrl_cell_ids <- colnames(ctrl)
cre_cell_ids_base <- cre_cell_ids %>% str_split(pattern = "-", simplify = TRUE) %>% .[ , 1]


DimPlot(object = ctrl,
        reduction = "umap",
        label = TRUE) + NoLegend()

ko <- obj.list[["ko"]]
ko_cell_ids <- colnames(ko)

DimPlot(object = ko,
        reduction = "umap",
        label = TRUE) + NoLegend()
```

3. Bring in loom files for ctrl and ko samples

```
# Use devtools to install hdf5r and loomR from GitHub
# devtools::install_github(repo = "hhoeflin/hdf5r")
# devtools::install_github(repo = "mojaveazure/loomR"
# library(devtools)
# install_github("velocyto-team/velocyto.R")   
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github("mojaveazure/seurat-disk")

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(loomR)

### Load filtered Seurat object and loom file
seurat_object <- readRDS('/Location/of/seurat_object.rds')
ldat <- ReadVelocity(file = "/Location/of/loom_file.loom")

for (i in names(x = ldat)) {
  ### Store assay in a new variable
  assay <- ldat[[i]]
  
  ### Rename cell names in loom file to match cell names in Seurat object
  
  colnames(assay)[str_detect(colnames(assay), pattern = "A1_")] <- colnames(assay)[str_detect(colnames(assay), pattern = "A1_")] %>% str_replace("x$", "-1_1")
  colnames(assay)[str_detect(colnames(assay), pattern = "A2_")] <- colnames(assay)[str_detect(colnames(assay), pattern = "A2_")] %>% str_replace("x$", "-1_2")
  colnames(assay) <- gsub('Part of cell name to change', 'Changed part', colnames(assay))
  colnames(assay) <- gsub('A1_CKDL210009739-1a-SI_TT_B3_HC2W5DSX2:', '', colnames(assay))
  colnames(assay) <- gsub('A2_CKDL210009740-1a-SI_TT_B6_HC2W5DSX2:', '', colnames(assay))
  
  ### Subset to filtered cells in Seurat object
  assay <- assay[,colnames(seurat_object)]
  
  ### Add assay to Seurat object
  seurat_object[[i]] <- CreateAssayObject(counts = assay)
}

DefaultAssay(seurat_object) <- "RNA"
SaveH5Seurat(seurat_object, filename = "cre.h5Seurat")
Convert("cre.h5Seurat", dest = "h5ad")

