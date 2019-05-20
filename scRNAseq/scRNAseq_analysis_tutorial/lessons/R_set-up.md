# Set-up for using R on O2 for single-cell RNA-seq

Log onto O2 and start an interactive session - more memory is often required, and I start with `--mem` of `64G` for ~ 10,000 - 15,000 cells. It also helps to have X11 forwarding enabled, `--x11`, to work through the analysis and view plots.

If X11 is not working for you, then make sure you follow the directions given on the [O2 wiki](https://wiki.rc.hms.harvard.edu/display/O2/Using+X11+Applications+Remotely).

Generally, I have my script open in a different terminal window and start an interactive session to copy and paste each line of the script. I explore the output and decide whether to adjust parameters. Alternatively, you could just run the script and explore the output.

```r
srun --pty -p interactive -t 0-12:00 --x11 --mem 64G /bin/bash
```

We should specify the R library to use in our `~/.Renviron` file. You can either create a personal R library or use those provided by the core. Either way you will need to store the path to the `R_LIBS_USER` variable. 

Within the same file, it will help to set the `R_MAX_NUM_DLLS` variable to a high number in order to use many of the single cell packages. I have mine set to `200`, which has worked so far.

```r
# Library using R 3.5.1 with Seurat 3.0
R_LIBS_USER="/n/data1/cores/bcbio/R/library/3.5.1-bioc-release_Seurat3.0"

# Library using R 3.5.1 but with version 2 Seurat
#R_LIBS_USER="/n/data1/cores/bcbio/R/library/3.5.1-bioc-release/library"

# Library using R 3.4
#R_LIBS_USER="/n/data1/cores/bcbio/R/library/3.4-bioc-release/library"

R_MAX_NUM_DLLS=200
```

Now to use R we need to load the required modules for our analysis:

```r
module load gcc/6.2.0 R/3.5.1 hdf5/1.10.1
```

> Be sure to match the module version of R to the library specified in the `~/.Renviron` file.
