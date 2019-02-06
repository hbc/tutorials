# Monocle

## Bringing in data

Bringing data into Monocle can be achieved in a few different ways depending on the source of the data. Regardless of source, a `CellDataSet` object needs to be created. A `CellDataSet` consists of:

- **expression matrix:** counts matrix, can be dense or sparse - in `CellDataSet` object stored in `assayData` slot
- **phenotype data:** metadata - in `CellDataSet` object stored in `pData` slot
- **feature data:** gene annotations, including a column named `gene_short_name` - in `CellDataSet` object stored in `fData` slot

Before we can construct the `CellDataSet`, the following libraries need to be loaded:

```r
library(Seurat)
library(monocle)
library(tidyverse)
library(readxl)
library(AnnotationHub)
library(ensembldb)
```

### Using Seurat output

Often pseudotime analysis is performed after performing QC, clustering and marker identification using Seurat. If this is the case, bringing the data into Monocle is quite easy. The `importCDS()` function is supposed to work, but I have had issues using this if my metadata is for a different number of cells than my raw data:

```r
# Read in Seurat object
seurat <- readRDS("path/to/seurat.rds")

# Create the 'CellDataSet'
cds <- importCDS(seurat)
```

Now we have our data stored as a `CellDataSet` object, and we can proceed through the Monocle workflow.

If this doesn't work for the creation of the object, then we can create the CDS using the metadata and raw counts.

### Using raw count matrix and metadata objects

**Step 1:** Read in count matrix - this should only be the filtered cells output from QC

If we have our Seurat object, we can access the different slots in which this data is stored:

```r
seurat <- readRDS("path/to/seurat.rds")

raw_counts <- seurat@raw.data

metadata <- seurat@meta.data
```

If no Seurat object, then we can read in the individual count matrix and metadata objects: 

```r
# Bring in count matrix from bcbio
raw_counts <- readMM("path/to/tagcounts.mtx")

# Assign row names and column names of matrix
gene_names <- read.csv("path/to/tagcounts.mtx.rownames", header = FALSE)

cell_ids <- read.csv("path/to/tagcounts.mtx.colnames", header = FALSE)

rownames(raw_counts) <- gene_names[, 1]

colnames(raw_counts) <- cell_ids[, 1]
```

**Step 2:** Create the annotations - this could be brought in from the previous QC or clustering analysis or created again:

```r
# Acquire the gene names for the Ensembl IDs
## Connect to AnnotationHub
ah <- AnnotationHub()

## Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

## Acquire the latest annotation files
id <- ahDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)

## Download the appropriate Ensembldb database
edb <- ah[[id]]

## Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

## Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
```

Now that we have the annotations, we can subset those to only the genes present in the counts data:

```r        
## Subset to include only those genes in data frame
monocle_annotations <- annotations[which(annotations$gene_id %in% rownames(raw_counts)), ]
```

Monocle expects the annotations to be formatted with the gene IDs as row names and a column corresponding to gene symbol named `gene_short_name`:

```r
## Make the row names of the annotations to be the same gene IDs as in the count matrix
rownames(monocle_annotations) <- monocle_annotations$gene_id

## Change name of gene symbol column to 'gene_short_name'
colnames(monocle_annotations)[colnames(monocle_annotations) == "gene_name"] <- "gene_short_name"
```

Finally, the order of genes need to match between the features and the counts:

```r
# Check if all genes are annotated
which(!(rownames(raw_counts) %in% rownames(monocle_annotations)))

# Remove genes not annotated
raw_counts <- raw_counts[which(rownames(raw_counts) %in% rownames(monocle_annotations)), ]

## Check all of the row names of the annotations match the row names of the counts
all(rownames(raw_counts) %in% rownames(monocle_annotations))

all(rownames(raw_counts) == rownames(monocle_annotations))

## If not, then match them
idx <- match(rownames(raw_counts), rownames(monocle_annotations))

monocle_annotations <- monocle_annotations[idx, ]

## Sanity check
all(rownames(raw_counts) == rownames(monocle_annotations))
```

Then, we can create the feature data object used to create the `CellDataSet` as an `AnnotatedDataFrame`:

```r
## Create feature data as an 'AnnotatedDataFrame'
fd <- new("AnnotatedDataFrame", data = monocle_annotations)
```

**Step 3:** Create an `AnnotatedDataFrame` to be used to create the `CellDataSet`:

```r
## Read in the metadata if not already present
metadata <- read.csv("path/to/metadata.csv")

## Check that the columns of the counts corresponds to the rows of the metadata
all(rownames(metadata) == colnames(raw_counts))

# if not matching, then use match() similar to above

## Create the phenotype data as an 'AnnotatedDataFrame'
pd <- new("AnnotatedDataFrame", data = metadata)
```

**Step 4:** Create the `CellDataSet` object - `expressionFamily` depends on the type of data. 

```r
cds <- newCellDataSet(raw_counts,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())
```

> **NOTE:** If the data have UMIs, and it's not an extremely small dataset, then `negbinomial.size()` is the correct option. More details available in the [Monocle docs](http://cole-trapnell-lab.github.io/monocle-release/docs/#choosing-a-distribution-for-your-data-required).

> **NOTE:** ## Using Cell Ranger output
>
>Taken directly from the Monocle documentation: 'If you have 10X Genomics data and are using cellrangerRkit, you can use it to load your data and then pass that to Monocle as follows:'
>
>```r
>cellranger_pipestance_path <- "/path/to/your/pipeline/output/directory"
>gbm <- load_cellranger_matrix(cellranger_pipestance_path)
>
>fd <- fData(gbm)
>
># The number 2 is picked arbitrarily in the line below.
># Where "2" is placed you should place the column number that corresponds to your
># featureData's gene short names.
>
>colnames(fd)[2] <- "gene_short_name"
>
>gbm_cds <- newCellDataSet(exprs(gbm),
>                  phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
>                  featureData = new("AnnotatedDataFrame", data = fd),
>                  lowerDetectionLimit = 0.5,
>                  expressionFamily = negbinomial.size())
>```

## Estimating size factors and dispersions

Similar to any other RNA-seq analysis exploring differential expression, we need to calculate the size factors for normalization and dispersions per gene:

```r
# Estimate the size factors
cds <- estimateSizeFactors(cds)

# Estimate the gene dispersions
cds <- estimateDispersions(cds)
```

## Additional QC suggested by Monocle

The Monocle tutorial suggests filtering low quality cells for minimum expression levels and for doublets. We likely have already performed the filtering for minimum expression during the original QC, but the doublet filtering has not been performed. Evidently, the trajectory analysis is quite sensitive to the presence of doublets. Moving forward we may consider different tools to perform this filtering.

```r
# Additional filtering of low quality cells

# Removing lowly expressed genes
cds <- detectGenes(cds, min_expr = 0.1)

# Removing genes not expressed in at least 10 cells
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# Removing cells that may be doublets
pData(cds)$Total_mRNAs <- Matrix::colSums(exprs(cds))

cds <- cds[,pData(cds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) +
            2*sd(log10(pData(cds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(cds)$Total_mRNAs)) -
            2*sd(log10(pData(cds)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(cds), color = viralLoad, geom =
"density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)

cds <- cds[,pData(cds)$Total_mRNAs > lower_bound &
      pData(cds)$Total_mRNAs < upper_bound]
cds <- detectGenes(cds, min_expr = 0.1)
```

After performing the filtering, it is suggested to check the data to ensure the counts follow an approximate log-normal scale.

```r
# Check expression to make sure filtered counts follow approximate log-normal distribution

# Log-transform each value in the expression matrix.
L <- log(exprs(cds[expressed_genes,]) + 1)

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

# Plot the distribution of the standardized gene expression values.
qplot(value, geom = "density", data = melted_dens_df) +
stat_function(fun = dnorm, size = 0.5, color = 'red') +
xlab("Standardized log(norm_counts)") +
ylab("Density")
```

Now we are ready for classifying cells by cell type. 

## Cell classification

Monocle uses a bit different method for clustering cells by taking in known marker genes to aid with clustering and identification of cell type. We need to provide Monocle with the gene IDs for the marker genes of the different clusters.

For example, if working with immune cells, we could have identified good cell markers for our dataset with Seurat previously:

```r
# Acquiring the rownames of markers
CD14_id <- row.names(subset(fData(cds), gene_short_name == "CD14")) # Monocytes
CD3_id <- row.names(subset(fData(cds),
                             gene_short_name == "CD3D")) # T cells
CD4_id <- rownames(subset(fData(cds),gene_short_name == "CD4")) # CD4+ T cells
CD8_id <- row.names(subset(fData(cds),
                           gene_short_name == "CD8A")) # CD8+ T cells
CD19_id <- row.names(subset(fData(cds),
                           gene_short_name == "CD19")) # B cells

# Creating hierarchy for cell assignment
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Monocytes", classify_func =
                           function(x) { x[CD14_id,] > 1 })
cth <- addCellType(cth, "T cell", 
                   classify_func=function(x) {x[CD3_id,] > 0})

cth <- addCellType(cth, "CD4+ T cells", classify_func = function(x)
{ x[CD14_id,] < 1  & x[CD4_id,] > 1 &  x[CD8_id,] < 1 }, 
parent_cell_type_name = "T cell")

cth <- addCellType(cth, "CD8+ T cells", classify_func = function(x)
{ x[CD14_id,] < 1  & x[CD8_id,] > 1 &  x[CD4_id,] < 1 }, 
parent_cell_type_name = "T cell")

cth <- addCellType(cth, "B cells", classify_func =
                           function(x) { x[CD19_id,] > 1 })

```

Now that we have the heirarchy for cell type assignment, we can assign cells to a known cell type. A cell is assigned to a cell type if at least the fraction of counts specified with the `frequency_thres` argument correspond to that cell type marker.

```r
# Should assign cells to one of the cell types specified in the heirarchy, Ambiguous, or Unknown 
cmv <- classifyCells(cds = cds, cth = cth, frequency_thres = 0.1)
``` 

We can explore the assignments to see if they make sense. At this stage in the analysis, it is normal for the majority of cells to be of 'Unknown' cell type. However, if there are a lot of 'Ambiguous' cells, then you may want to modify your assignment heirarchy.

```r
# Check number of cells per celltype - at this stage majority of cells are often unknown
table(pData(cmv)$CellType)

# Visualize by pie chart
pie <- ggplot(pData(cmv),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)

pie + coord_polar(theta = "y") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Check for specific cell types if desired
subset(pData(cmv), CellType == "pDCs")
```

## Identify clustering genes using an 'Unsupervised' method

Now we can try to assign identity to the 'Unknown' cells by using the prinicipal components that explain the largest amount of variance in the data, somewhat similar to Seurat's method.

```r
# Assign celltype to Unknown cells
disp_table <- dispersionTable(cmv)

# Identify ordering genes - unsupervised clustering

## Subset those genes with expression higher than 0.1
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

## Mark genes to be used for clustering
cmv <- setOrderingFilter(cmv, unsup_clustering_genes$gene_id)

## View genes to be used for clustering
plot_ordering_genes(cmv)

## Determine number of principal components to use based on where elbow meets the surface
# x11(type="cairo") # Run if error viewing the following plot
plot_pc_variance_explained(cmv, return_all = F) # norm_method='log'
```

Now we can perform the dimensionality reduction using the identified principal components and the tSNE method. You will need to choose the number of clusters to return; I randomly chose 15 clusters, but you could choose more or less based on expectations. If you choose more, than expected, you can always merge together later on in the analysis. Also, we need to choose the number of dimensions to use, which should be based on analysis of the Elbow (skree) plot where the elbow just seems to touch the base.

```r
## Reduce dimensions for tSNE viewing with max components of 2 and number of dimensions equal to the PCs determined in elbow plot
cmv <- reduceDimension(cmv, max_components = 2, num_dim = 9,
                reduction_method = 'tSNE', verbose = T)

## Cluster the cells to a certain number of clusters - will limit # clusters returned - randomly chose 15, but may return less                
cmv <- clusterCells(cmv, num_clusters = 15)
```

Let's explore the quality of our clustering by checking our known markers:

```r
## Explore cluster assignment
head(pData(cmv))

plot_cell_clusters(cmv, 1, 2, color = "CellType",
    markers = c("CD14", "CD36", "CD3D", "CD8A", "CD4", "CD19"))
    
cmv <- reduceDimension(cmv, max_components = 2, num_dim = 9,
            reduction_method = 'tSNE',
            residualModelFormulaStr = "~ condition + num_genes_expressed",
            verbose = T)

cmv <- clusterCells(cmv, num_clusters = 15)

plot_cell_clusters(cmv, 1, 2, color = "Cluster") +
    facet_wrap(~CellType)
 ```
 
## Further identify clustering genes using a 'Supervised' method

While the unsupervised clustering method allowed for using genes that were more highly expressed and variable for determining the principal components to use for clustering, the supervised method will instead choose genes that co-vary with the cell type markers given. After identifying the genes that co-vary signficantly with the cell type markers, we will select genes with high specificity; usually it's best to pick the top 10 or 20 genes most specific per cell type.

```r
# Identify ordering genes - supervised method

# Identifying genes that co-vary with markers
marker_diff <- markerDiffTable(cmv[expressed_genes,],
            cth,
            residualModelFormulaStr = "~ condition + num_genes_expressed",
            cores = 1)

# Selecting the genes that significantly co-vary
candidate_clustering_genes <-
    row.names(subset(marker_diff, qval < 0.01))

# Determine specificity of the markers
marker_spec <- calculateMarkerSpecificity(cmv[candidate_clustering_genes,], cth)

head(selectTopMarkers(marker_spec, 3))
```

Now we can use these specific markers to cluster the cells. We will pick the top 500 markers for each cell type (although I am unsure why the choice is 500 here and not the 10-20 genes mentioned previously. We determine the unique top 500 specific markers for each cluster, then mark the genes 
```r
# Select the specific cell type genes to use
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)

# Mark that these are the genes to be used for clustering 
cmv <- setOrderingFilter(cmv, semisup_clustering_genes)

# Explore the variance explained by the genes
plot_ordering_genes(cmv)

plot_pc_variance_explained(cmv, return_all = F)

# Use these genes for the clustering
cmv <- reduceDimension(cmv, max_components = 2, num_dim = 9,
  norm_method = 'log',
  reduction_method = 'tSNE',
  residualModelFormulaStr = "~ condition + num_genes_expressed",
  verbose = T)
 
# Cluster the genes similar to previously
cmv <- clusterCells(cmv, num_clusters = 15)

# Explore the clustering
plot_cell_clusters(cmv, 1, 2, color = "CellType",
    markers = c("CD14", "CD36", "CD3D", "CD8A", "CD4", "CD19"))
    
plot_cell_clusters(cmv, 1, 2, color = "Cluster") +
    facet_wrap(~CellType)
```

## Imputing cell types

For those cells that are still of 'Unknown' cell type, we can impute the identity based on the expression of markers from the other cells in that cluster. We will impute the identities of the 'Unknown' cells using a threshold of 10% for the percentage of cluster marked as a certain type of cell to impute the values of the remaining cells.

```r
# Impute cell type
imputed <- clusterCells(cmv,
              num_clusters = 15,
              frequency_thresh = 0.1,
              cell_type_hierarchy = cth)
              
plot_cell_clusters(imputed, 1, 2, color = "CellType",
    markers = c("CD14", "CD36", "CD3D", "CD8A", "CD4", "CD19"))
    
pie <- ggplot(pData(imputed),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)

pie + coord_polar(theta = "y") +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())

table(pData(imputed)$CellType)
```

## Trajectory analysis

To perform trajectory analysis, you will need to subset out the cells of interest from your object. In this example, I'm interested in monocytes.

```r
# Subset out monocytes
monocyte_cells <- row.names(subset(pData(imputed), CellType == "Monocytes"))

monocytes <- imputed[ , monocyte_cells]

dim(monocytes)
```

```r
# Trajectory analysis
monocytes <- detectGenes(monocytes, min_expr = 0.1)

```
