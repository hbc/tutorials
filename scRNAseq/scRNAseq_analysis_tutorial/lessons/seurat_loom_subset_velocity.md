# Subset loom file to cells of interest

```r
# Get loom data
ldat <- ReadVelocity(file = "data/merged_loom_files/all_merged.loom")

# Turn loom to seurat
bm <- as.Seurat(x = ldat)
bm[["RNA"]] <- bm[["spliced"]]

# Bring in Seurat with known clusters
seurat_object <- readRDS('data/seurat_combined_sct_umap_FBS.rds')

# Change cell ids of loom seurat to match cell ids in seurat with known clusters
all_cells <- Cells(bm)
all_cells[str_detect(all_cells, pattern = "A1_")] <- all_cells[str_detect(all_cells, pattern = "A1_")] %>% str_replace("x$", "-1_1")
all_cells[str_detect(all_cells, pattern = "A2_")] <- all_cells[str_detect(all_cells, pattern = "A2_")] %>% str_replace("x$", "-1_2")
all_cells[str_detect(all_cells, pattern = "A3_")] <- all_cells[str_detect(all_cells, pattern = "A3_")] %>% str_replace("x$", "-1_3")
all_cells[str_detect(all_cells, pattern = "A4_")] <- all_cells[str_detect(all_cells, pattern = "A4_")] %>% str_replace("x$", "-1_4")

all_cells <- gsub('A1_CKDL210009739-1a-SI_TT_B3_HC2W5DSX2:', '', all_cells)
all_cells <- gsub('A2_CKDL210009740-1a-SI_TT_B6_HC2W5DSX2:', '', all_cells)
all_cells <- gsub('A3_CKDL210009741-1a-SI_TT_B2_HC2W5DSX2:', '', all_cells)
all_cells <- gsub('A4_CKDL210009742-1a-SI_TT_B7_HC2W5DSX2:', '', all_cells)

new_names <- all_cells
bm <- RenameCells(bm, new.names = new_names)

# Get names of seurat_object
DefaultAssay(seurat_object) <- "RNA"

sub_genes <- rownames(seurat_object)
sub_cells <- colnames(seurat_object)

# Subset loom seurat
bm <- subset(bm, features = sub_genes, cells = sub_cells)

# Add  cluster ID to metadata file for each cell
bm <- AddMetaData(bm, seurat_object@meta.data[, c("DE_group", "sample_simple", "seurat_clusters")])

# Add all slots to object
bm@reductions[["pca"]] <- seurat_object@reductions[["pca"]]
bm@reductions[["umap"]] <- seurat_object@reductions[["umap"]]

bm@assays$integrated <- seurat_object@assays$integrated
bm@assays$SCT <- seurat_object@assays$SCT


# Split into individual and merged objects

# Save object and convert to h5ad format for scvelo

DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "all_samples.h5Seurat")
Convert("all_samples.h5Seurat", dest = "h5ad")

```

#################

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

