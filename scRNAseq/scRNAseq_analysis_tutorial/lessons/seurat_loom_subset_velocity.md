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

library(loomR)

crtl_sub_loom <- subset(ctrl_loom_merged, m = ctrl_cell_ids, n = NULL, filename = NULL,
  chunk.size = 1000, overwrite = FALSE, display.progress = TRUE)
  
ldat <- ReadVelocity(file = "path_to_ctrl.loom")
bm <- as.Seurat(x = ldat)

