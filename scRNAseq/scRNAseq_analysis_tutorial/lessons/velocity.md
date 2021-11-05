# Velocity analysis using velocyto, Seurat and scVelo

> _**NOTE:** To run analysis in a Jupyter notebook on O2, please follow the [available instructions](https://github.com/hbc/knowledgebase/blob/master/rc/jupyter_notebooks.md)._

1. Generate loom files containing layers for spliced and unspliced reads using Velocyto.

  ```
  # On command line
  
  # Load modules
  module load gcc/9.2.0 python/3.8.12
  
  # Create virtual environment
  # virtualenv velocyto --system-site-packages
  
  # Activate virtual environment
  source velocyto/bin/activate

  # Install tools
  # pip3 install numpy scipy cython numba matplotlib scikit-learn h5py click
  # pip3 install velocyto
  # pip3 install scvelo
  
  # Run velocyto
  velocyto run10x -m ../data/mm10_rmsk.gtf ../final/cellranger_6.0.0/count/expect_cells/A1_CKDL210009739-1a-SI_TT_B3_HC2W5DSX2/ genes.gtf
  
  # Create merged loom object - start by copying first file in 'files' below, then add another loom file in python below.
  # cp path_to_file1.loom path_to_file1_2_merged.loom
  ```

2. In python merge any loom files desired

  ```python
  #python3
  import velocyto as vcy
  import loompy
  import scvelo as scv
  import numpy as np
  import h5py
  import scipy
  import cython
  import numba
  import matplotlib
  import 
  import click
  
  # files = ["path_to_file1.loom", "path_to_file2.loom"]
  files = ["../final/cellranger_6.0.0/count/expect_cells/A1_CKDL210009739-1a-SI_TT_B3_HC2W5DSX2/velocyto/A1_CKDL210009739-1a-SI_TT_B3_HC2W5DSX2.loom", "../final/cellranger_6.0.0/count/expect_cells/A2_CKDL210009740-1a-SI_TT_B6_HC2W5DSX2/velocyto/A2_CKDL210009740-1a-SI_TT_B6_HC2W5DSX2.loom", "../final/cellranger_6.0.0/count/expect_cells/A3_CKDL210009741-1a-SI_TT_B2_HC2W5DSX2/velocyto/A3_CKDL210009741-1a-SI_TT_B2_HC2W5DSX2.loom", "../final/cellranger_6.0.0/count/expect_cells/A4_CKDL210009742-1a-SI_TT_B7_HC2W5DSX2/velocyto/A4_CKDL210009742-1a-SI_TT_B7_HC2W5DSX2.loom"]

  ds = loompy.connect("data/merged_loom_files/all_merged.loom")
  for fn in files[1:]:
    ds.add_loom(fn, batch_size=1000)
  ```
  
3. Perform all QC, normalization and clustering using Seurat as described at [http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html). 

    If you cannot do your own analyses and must use the data given to you by the client, then documentation for merging/subsetting loom files with their Seurat data/clusters is [available](https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/seurat_loom_subset_velocity.md).
  
    ```r
    deactivate velocyto
    module load gcc/6.2.0 R/4.1.1
    
    # library(devtools)
    # install_github("velocyto-team/velocyto.R")   
    # remotes::install_github('satijalab/seurat-wrappers')
    # remotes::install_github("mojaveazure/seurat-disk")
    
    library(Seurat)
    library(SeuratDisk)
    library(SeuratWrappers)
    
    ldat <- ReadVelocity(file = "path_to_file.loom")
    bm <- as.Seurat(x = ldat)
    bm[["RNA"]] <- bm[["spliced"]]
    bm <- SCTransform(bm)
    bm <- RunPCA(bm)
    bm <- RunUMAP(bm, dims = 1:30)
    bm <- FindNeighbors(bm, dims = 1:30)
    bm <- FindClusters(bm)
    DefaultAssay(bm) <- "RNA"
    SaveH5Seurat(bm, filename = "mouseBM.h5Seurat")
    Convert("mouseBM.h5Seurat", dest = "h5ad")
    ```
  
4. Use scVelo in python to construct velocity estimates and trajectories and continue following [http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html](http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html) and consulting [scVelo documentation](https://scvelo.readthedocs.io/VelocityBasics)

  ```python
  import scvelo as scv
  adata = scv.read("mouseBM.h5ad")
  adata
  
  scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=3000)
  
  scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
  
  scv.tl.velocity(adata)
  
  scv.tl.velocity_graph(adata)
  
  scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")
  
  scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120)
  
  scv.tl.recover_dynamics(adata)
  
  scv.tl.latent_time(adata)
  
  scv.pl.scatter(adata, color="latent_time", color_map="gnuplot")
  
  top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
  scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time", col_color="seurat_clusters", n_convolve=100)
  ```
