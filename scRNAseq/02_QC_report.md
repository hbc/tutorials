# bcbioSingleCell QC Report

### Different approaches to running `bcbioSingleCell`

There are various approaches to running `bcbioSingleCell` to generate the QC report. The first part is getting all of the output from the **`bcbio` final directory loaded in to create the `bcb` object**. The next step is running through code which will **compute metrics and generate figures for quality assessment**. This second step is best done locally so you can run the code interactively and assess things as you run through the report code chunks. For approaches #1 and #2 listed below, you are doing everything locally. For #3 and #4 you are creating the `bcb` object on the cluster, and then moving it local to create the report. 

1. **Using a [Docker image](https://hub.docker.com/r/lpantano/bcbiosinglecell/)**. 
    - First, install Docker 
    - Pull the Docker image: `docker pull lpantano/bcbiosinglecell:r3.5-bsc0.1.5`
    - Set your memory RAM limit to 4G or more. This is done with the Docker main application (Preferences -> Advanced)
    - Mount the O2 `final` directory from your `bcbio` run on your laptop
    - Open up a terminal and make sure you are in your home directory (or a place where you can easily navigate to the mount space)
    - Run the Docker image: `docker run -d -p 8787:8787 -e ROOT=TRUE -v $(pwd):/home/rstudio lpantano/bcbiosinglecell`
    - In a browser connect to RStudio: localhost:8787 with user and password: rstudio/rstudio. From here you can start [Creating the metadata file](#metadata), and continue working within the Docker container to create the QC report.   
    
    ---
    
    > **NOTE:** If you start a Docker container and realize you want to start a new one, you will want to kill this one and remove it using the commands below:
    > ```
    > docker ps # to see your containers listed by id
    > docker stop <container_id>
    > docker rm <container_id> ```

2. Running it **locally on your laptop RStudio**. This will require you to install `bcbioSingleCell` and also mount O2. Note that if you have more 400K-500K cells this will max out of memory. Also, note you may have to deal with problems with various dependency packages as you update R. If you choose this method, skip down to [Creating the metadata file section](#metadata) and get started.

3. **Generate the `bcb` object on the O2 cluster**. You are limited to using R 3.4.1 because that is what is available for conda and the modules, but `bcbioSingleCell` is backwards compatible to R 3.4.1. The code is as follows:

```r
	bcbio <- loadSingleCell("~/bcbio/PIs/path/to/final/",
                        interestingGroups = "sampleName",
                        sampleMetadataFile = "~/path/to/metadata", 
                        gtfFile = "~/bcbio/PIs/path/to/Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf")
	
	save(bcbio_output, file="data/bcb.rda")
```

The above code chunk can be run on O2 in one of two ways:

   - **A.** Using a **conda install of R 3.4.1** and pointing to a [shared R library](#rlib). For the conda recipe you can find more information [here](https://steinbaugh.com/r_bioconda). Keep note of the different versions when you create your environment (i.e. pandoc 1 is required for rmarkdown (version 2 is super buggy) and hdf5 1.10.1 is required for the latest version of Seurat, or it won’t compile)
	
   - **B.** Using the **R 3.4.1 module** and pointing to the [shared R library](#rlib). This may require some troubleshooting with the HMSRC folks as it has been known to be problematic.
	
		
> #### Using a pre-existing shared R library on O2 (for single cell RNA-seq) <a name="rlib"></a>
>  This library has been created for use with single cell RNA-seq analysis. It can be used not only for QC but also for clustering with Seurat. First, you will need to edit your `.Renviron` file to have the following inside:
> 
> ```
> R_LIBS_USER="/n/data1/cores/bcbio/R/library/3.4-bioc-release/library"
> R_MAX_NUM_DLLS=150
> ```
> 
> Then start an interactive session with extra memory and x11:
> 
> `$ srun --pty -p interactive -t 0-12:00 --x11 --mem 128G /bin/bash`
> 
> After starting the interactive session, load the necessary R modules and start R as described at https://github.com/hbc/knowledgebase/blob/master/research/scrnaseq/Single-Cell.md#shared-installation-in-o2.



### Creating the metadata file <a name="metadata"></a>

Use the information from the client to construct the metadata table to use with bcbioSingleCell R package according to the specifications detailed at [https://github.com/hbc/bcbioSingleCell](https://github.com/hbc/bcbioSingleCell). You will need the columns for `description`, `index`, `sequence`, and `sampleName`. You can add any additional metadata as desired.

- **Example metadata table:**
	
	![example metadata](../img/sc_metadata.png)
	
	- **Important:** the `sequence` column for the inDrop metadata is the **Forward** sequence, not the same as the sequences present in the `sample_barcodes` file, which is the reverse complement. 

### Quality control report <a name="qc"></a>

#### Setting up

1. Choose the quality control template.

	> Documentation for all functions available from the bcbioSingleCell package is available at [http://bioinformatics.sph.harvard.edu/bcbioSingleCell/reference/index.html](http://bioinformatics.sph.harvard.edu/bcbioSingleCell/reference/index.html)

2. Edit the information in the files `_header.Rmd` and `_footer.Rmd` with experiment-specific information.

3. Install `bcbioSingleCell` and load the library:
	
	```r
	# devtools::install_github("hbc/bcbioSingleCell") # Add argument `ref = "develop"` if need development branch
	
	library(bcbioSingleCell)
	```
	
4. Bring in data from bcbio:
	
	```r
	bcbio <- bcbioSingleCell("~/bcbio/PIs/path/to/final/",
                    organism = "Homo sapiens",
                    interestingGroups = "sampleName",
                    sampleMetadataFile = "~/path/to/metadata.csv",
                    ensemblRelease = 92L,
                    genomeBuild = "GRCh38")
	
	save(bcbio_output, file="data/bcb.rda")
	```


5. Choose the filtering parameters to use. You can start with these parameters, then after viewing the data, change to better values. Generally, you don't want `minGenes`/`minUMIs` to be any lower than 500.  You would hope for at least 1000 genes/UMIs detected per sample. After choosing parameters, run the entire `r setup` chunk by clicking on the green triangle at the top of the setup chunk (if you clear your environment, you need to run the chunk this way to make the `params` reappear.
	
	**Choosing parameters**
	```r
	params:
    	  bcb_file: "data/bcb.rda"
    	  min_genes: 500
          max_genes: !r Inf
          max_mito_ratio: 0.25
          min_novelty: 0.85
          min_cells_per_gene: 10
          data_dir: !r file.path("data", Sys.Date())
  	```
	
	**Running setup chunk**
	```r
	# Shared RMarkdown settings
	prepareSingleCellTemplate()
	if (file.exists("setup.R")) {
	    source("setup.R")
	}

	# Directory paths
	dataDir <- file.path(params$outputDir, "data")

	# Load bcbioSingleCell object
	bcbName <- load(params$bcbFile)
	bcb <- get(bcbName, inherits = FALSE)
	```
	
	```r
	eval=file.exists("_header.Rmd")
	```

	```r
	sampleMetadata(bcb)
	```

6. For the count alignment, be sure to update the **linked Ensembl** to be accurate for the organism. This information is present in the file: `_footer.Rmd`. 

7. To explore the raw data stored inside the `bcb` object, the following functions can be helpful:
	
	```r
	# Access metadata for each sample: "sampleID", "sampleName", "description", "fileName", "index", "sequence", "revcomp"
	sampleMetadata(bcb)
	
	# Access metadata for each cell: "nCount", "nUMI", "nGene", "nCoding", "nMito", "log10GenesPerUMI", "mitoRatio" 
	colData(bcb)
	
	# Access raw counts - each column represents a single cell
	counts <- counts(bcb)
	
	# Can return cells from a particular sample by using metadata information about which sample corresponds to each barcode
	unsort_counts <- counts[, str_detect(colnames(counts), "run1_ATTAGACG")] # Return only the counts for the `Unsorted` sample
	
	# Extract information associated with each gene including "ensgene", "symbol", "description", "biotype", "broadClass"
	rowData(bcb)
	
	# Return the genes that are associated with a broad class (ex: mitochondrial contamination)
	subset(rowData(bcb), broadClass == "mito")
	```

#### Quality Control Metrics

##### Reads per cell

8. Evaluate the number of reads per cell:

	```r
	plotReadsPerCell(bcb)
	```
	
	The three plots give different ways of looking at the number of reads per cell. Generally you would like to see a large peak at around 10,000 reads per cell, and you hope your filtering threshold of 1,000 reads per cell used in bcbio has removed the poor quality cells with few number of reads. The filtering threshold of 1,000 is represented by the vertical dotted line.

	For example, in the figures below, the yellow sample is worrisome because we see a small peak at 10,000 reads per cell, but a much larger peak at 1,000 reads per cell. The larger peak merges into the poor quality cells with few reads per cell.
	
	<img src="../img/sc_qc_reads_ridgeline.png" width="600">
	
	The proportional histogram looks a bit better, as you hope to see all of the samples with peaks in relatively the same location between 10,000 and 100,000 reads per cell. However, the yellow sample still has this shoulder, which is indicative of many poor quality cells. If this were the only issue with the data, we may want to set the threshold to be more strict to ~10,000 reads per cell to get rid of the cells constituting the shoulder in the yellow sample.

	<img src="../img/sc_qc_reads_histogram.png" width="500">
	
##### Cell counts

9. Determine the number of cells detected per sample:

	```r
	plotCellCounts(bcb)
	```

	The cell counts are determined by the number of unique cellular barcodes detected. During the inDrop protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets with a single cell and lysis/reaction mixture. Upon treatment of UV and cell lysis, all components mix together inside the droplet and reverse transcription proceeds, followed by droplet breakup and linear amplification for library preparation. While each hydrogel should have a single cellular barcode associated with it, occasionally a hydrogel can have more than one cellular barcode. We often see all possible combinations of cellular barcodes at a low level, leading to a higher number of cellular barcodes than cells.

	You expect the number of unique cellular barcodes to be around the number of sequenced cells (determined in step 1) or greater due to some hydrogels having more than one cellular barcode. The yellow sample below seems to have at least double the number of cellular barcodes as the other samples.

	<img src="../img/sc_qc_cellcounts.png" width="500">

##### UMI counts per cell

10. Determine the number of UMI counts (transcripts) per cell:

	```r
	plotUMIsPerCell(
    		bcb,
	    	min = params$minUMIs)
	```

	The UMI counts per cell should be generally above 500, although usable, it's still low if between 500-1000 counts. If UMIs per cell is 500-1000 counts, then the cells probably should have been sequenced more deeply. The threshold of 500 was given in the `params`, and this is represented by the vertical dashed line in the plots.
	
	The number of UMIs per cell tends to be very low for the Unsorted sample (yellow). The other samples have good numbers of UMIs per cell, indicating a problem only with the Unsorted sample. Using this cutoff, we will lose the majority of the Unsorted cells.
	
	<img src="../img/sc_qc_umisPerCell.png" width="500">
	
##### Genes detected per cell

11. Discover the number of genes detected per cell:

	```r
	plotGenesPerCell(
	    bcb,
	    min = params$minGenes,
	    max = params$maxGenes)
	```

	Seeing gene detection in the range of 500-5000 is normal for inDrop analysis. Similar expectations for gene detection as for UMI detection.

	All samples other than the Unsorted sample have a good number of genes detected (with medians between 1,000 - 3,000 genes), which correspond to the numbers of UMIs per cell for each sample. However, the Unsorted sample has a very low median number of genes per cell, indicating a sample failure.

	<img src="../img/sc_qc_genesDetected.png" width="500">
	
##### UMIs vs. genes detected

12. Identify whether large number of poor quality cells present in any samples with low UMI/genes detected:

	```r
	plotUMIsVsGenes(bcb)
	```

	Poor quality cells are likely to have low genes and UMIs per cell. Therefore, a poor sample is likely to have cells in the lower left of the graph. Good cells should exhibit both higher number of genes per cell and higher numbers of UMIs. We also expect similar lines with similar slopes for all samples.
	
	The Unsorted sample has many cells with few UMIs and low number of genes per cell. The other samples look fine.
	
	<img src="../img/sc_qc_UMIsVsGenesDetected.png" width="500">
	
##### Mitochondrial counts ratio

13. Identify whether there is a large amount of mitochondrial contamination from dead or dying cells:

	```r
	plotMitoRatio(
	    bcb,
	    max = params$maxMitoRatio)
	```
	
	Poor quality samples for mitochondrial counts would have larger peaks above the 0.1 mitochondrial ratio mark, unless it is expected based on sample type. 
	
	There was just a very low number of genes detected for the Unsorted sample, so mitochondrial expression appears higher mainly due to this fact. The poor quality of the Unsorted sample does not appear to be due to dead or dying cells. The other samples have little mitochondrial expression, although hPSC sample has a bit more than the Sorted samples, and these cells will likely be removed using the threshold of 0.1. The hPSC sample was expected to contain brown adipocytes, which have higher quantities of mitochondrial expression, so it may have been advisable to keep these cells and move the threshold to 0.2.

	<img src="../img/sc_qc_mitoRatio.png" width="500">
	
##### Novelty

14. Explore the novelty for contamination with low complexity cell types:

	```r
	plotNovelty(
	    bcb,
	    min = params$minNovelty)
	```
	
	We can see the samples where we sequenced each cell less have a higher overall novelty, that is because we have not started saturated the sequencing for any given gene for these samples. Outlier cells in these samples might be cells that we have a less complex RNA species than other cells. Sometimes we can detect contamination with low complexity cell types like red blood cells via this metric.
	
	All of the samples look fine for complexity, except for the Unsorted sample, so it is unlikely that there is contamination with low complexity cell types in these of the samples. The Unsorted sample has a larger shoulder than desired, but is not bad by this metric.
	
	<img src="../img/sc_qc_novelty.png" width="500">
	

##### Filtered results

15. Run the filtering criteria and explore the plots again. The metrics should have improved greatly after removing low gene/UMI cells and high mitochondrial cells.

	```r
	bcbFiltered <- filterCells(bcb,
	minUMIs = params$minUMIs,
	minGenes = params$minGenes,
	maxGenes = params$maxGenes,
	maxMitoRatio = params$maxMitoRatio,
	minNovelty = params$minNovelty,
	minCellsPerGene = params$minCellsPerGene)
	```

	One main plot to look at to determine the success of the filtering criteria is the number of cell counts. You should expect roughly the number of sequenced cells per sample. We found out from the client that they had sequenced 2000-3000 cells, so the final numbers were around our expectations. If the number of cells sequenced is vastly different than the number returned after filtering, then you may need to re-visit the threshold criteria used for filtering.
	
	**Cell counts**
	
	<img src="../img/sc_qc_filtered_cellcounts.png" width="500">
	
	In addition, it is a good idea to explore all of the quality plots for the filtered data. All plots should be much improved for the number of reads per cell, genes detected, UMIs per cell, mitochondrial ratio, and novelty. The plots below show the filtered plots from the example data. Since the `Unsorted` sample was a poor quality sample, the filter will remove a large number of the cells for this sample; in this case all cells except 1 were filtered out. 
	
	**Reads per cell**
	
	The majority of cells have between 10,000 and 100,000 reads per cell, which is good.
	
	<img src="../img/sc_qc_filtered_reads.png" width="500">
	
	**Genes detected**
	
	The number of genes detected has also improved after the removal of the cells with low genes and or low UMIs.
	
	<img src="../img/sc_qc_filtered_genesDetected.png" width="500">
	
	**UMIs per cell**
	
	The numbers of UMIs per cell has also improved significantly, with the low quality cells dropped. It is worth noting here that the sample `Sort1` has many more UMIs per cell than the replicate `Sort2`. We will definitely want to regress out the variation due to numbers of UMI per cell in the clustering analysis.
	
	<img src="../img/sc_qc_filtered_umisPerCell.png" width="500">
	
	**UMIs versus genes detected**
	
	The correlations look more similar between samples, with few low gene and/or low UMI cells.
	
	<img src="../img/sc_qc_filtered_UMIsVsGenesDetected.png" width="500">
	
	**Mitochondrial ratio**
	
	The mitochondrial ratios are improved with no cells present with the high mitochondrial contamination.
	
	<img src="../img/sc_qc_filtered_mitoRatio.png" width="500">
	
	**Novelty**
	
	The novelty is also improved, with no shoulder for any of the samples.
	
	<img src="../img/sc_qc_filtered_novelty.png" width="500">
	
16. When you are satisfied with the filtered results, save the filtered data. You may need to adjust the filtering criteria multiple times to optimize the filtering results prior to saving the report.

	```r
	assignAndSaveData(name = "bcbFiltered", object = bcbFiltered, dir = dataDir)
	```
