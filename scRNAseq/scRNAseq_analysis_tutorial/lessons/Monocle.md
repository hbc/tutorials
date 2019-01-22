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

If no Seurat object, then we can read in the individual count matrix and metadata objects, which is a bit more complicated. 

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

Now we are ready for classifying cells by cell type. 

## Cell classification

Monocle uses a bit different method for clustering cells by taking in known marker genes to aid with clustering and identification of cell type.

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

# Unsupervised cell clustering

To perform clustering of the cells, we need to obtain the gene IDs for the genes that are expressed and the dispersions for these genes.

Then, we can identify those genes that have higher expression and dispersion values to use in the initial round of clustering.

```r
# Create column in fData with the number of cells expressed per gene
cmv <- detectGenes(cmv, min_expr = 1)

# Get the gene IDs for the genes that are expressed in more than 10 cells
expressed_genes <- row.names(subset(fData(cmv),
                                    num_cells_expressed >= 10))

# Extract the gene-wise dispersion values                                   
disp_table <- dispersionTable(cmv) 

# Return genes to use for unsupervised clustering by subsetting to genes with higher expression and dispersion
unsup_clustering_genes <- subset(disp_table,
                                 mean_expression >= 0.1) 
                                 
# Create column to denote whether to use gene to cluster cells
cmv <- setOrderingFilter(cmv, 
                         unsup_clustering_genes$gene_id)                                                                   
```
