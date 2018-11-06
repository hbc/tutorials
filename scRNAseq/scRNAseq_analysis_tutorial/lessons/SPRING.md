# SPRING

[SPRING](https://github.com/AllonKleinLab/SPRING) is a tool for visualizing and interacting with high dimensional data. The SPRING tool can be accessed [online](https://kleintools.hms.harvard.edu/tools/spring.html) or through a [local instance](https://github.com/AllonKleinLab/SPRING).

## Online webserver

When using the webserver for SPRING, the maximum number of cells to be used as input is **10,000 cells**. If looking to use SPRING for analysis of >10,000 cells, then a local instance would be required.

## Generating data for upload

Generally when using SPRING to visualize our data, we will need several files:

- the filtered raw counts
- gene symbols for genes in raw counts matrix
- metadata for the counts

### Filtered raw counts

The first file to attain is the filtered raw counts object. The easiest method is to access it from the Seurat object with clusters assigned. It is helpful if the clusters are named with the known or hypothesized cell types. So let's read in the seurat object from R:

```r
# Read in seurat object
seurat <- readRDS("path/to/data/seurat_tsne.rds")
```

**We can directly use this object if it contains less than 10,000 cells.** Alternatively, we could subset the object to the cells that we would like to view:

```r
# Get cell ids
spring_cells <- rownames(seurat@meta.data[seurat@meta.data$sample == "sample1" | seurat@meta.data$sample == "sample2", ])

# Subset the raw counts to these cells
spring_counts <- as.matrix(seurat@raw.data[ ,spring_cells])

# Write counts to file
write.csv(spring_counts, "spring/spring_counts.csv", quote= F)
```

> **NOTE:** We could subset by any other factor in the metadata similarly. We could also subset randomly to 10,000 cells by taking a sample of cells:
> 
> ```r
> # Get cell ids
> sampled_cells <- sample(x = seurat@cell.names, size = 10000, replace = F)
> 
> # Use cell ids to subset seurat
> spring_counts_10000 <- as.matrix(seurat@raw.data[ ,sampled_cells])
> ```

### Gene names for extracted counts

The next data we need is the gene names for the extracted counts. We can get the gene names directly from the counts file we just created.

```r
# Write genes to file
write(rownames(spring_counts), "spring/spring_genes.txt")
```


### Metadata

Finally, the last object we need is any metadata we might want to visualize. Now we can subset the Seurat object and extract the metadata stored in the `meta.data` slot of the Seurat object:

```r
# Add cell type to the metadata for each cell
seurat@meta.data$ident <- seurat@ident

# Use cell ids to subset seurat object
spring_seurat <- SubsetData(seurat, cells.use = spring_cells)

# Extract metadata including cell type
spring_meta <- spring_seurat@meta.data
```

#### Writing metadata to file

To write the metadata to file it needs to be in a particular format, which we can output using the `write()` function. We can specify any column of metadata that we would like to include in the SPRING visualizations.

```r
write(c("Cluster", spring_meta$ident), 
      file = "spring/spring_meta.csv", 
      sep = ",", 
      ncolumns = length(spring_meta$ident) + 1, 
      append = FALSE)

write(c("Condition", spring_meta$interestingGroups), 
      file = "spring/spring_meta.csv", 
      sep = ",", 
      ncolumns = length(spring_meta$ident) + 1, 
      append = TRUE)
      
write(c("Phase", spring_meta$phase), 
      file = "spring/spring_meta.csv", 
      sep = ",", 
      ncolumns = length(spring_meta$ident) + 1, 
      append = TRUE)
```

## SPRING interface

Now that we have the data that we would like to visualize, we can upload it to the [SPRING webserver](https://kleintools.hms.harvard.edu/tools/spring.html). 

Create a name for your dataset and a password, and choose the `Load new files` option. Load the following files by clicking on `Choose file`:

- **Expression data:** `spring_counts.csv`
- **Gene list:** `spring_genes.txt`
- **Cell groupings:** `spring_meta.csv`

Once loaded, select the `Upload` button. After the upload is successful, click on `Step 2: Process data`. There are a few parameters here that are available if you wish to adjust, then continue by clicking on `Begin processing data!`. Finally click on `Step 3: Click here to view data`, which should take you to the web interface for your data.

We encourage watching the 1-2 minute videos available on the [SPRING website](https://kleintools.hms.harvard.edu/tools/spring.html). Briefly, you can look at markers for different clusters by positively and negatively selecting cells. 

- Positive selection: `Shift`
- Negative selection: `Shift` + `Esc` 
- Deselection: `command`

After positively and negatively selecting cells, you can view the DE genes/marker genes by clicking on the `enriched genes` on the left-hand side.
