#!/bin/bash

#SBATCH -p priority             # partition name
#SBATCH -t 0-12:00              # hours:minutes runlimit after which job will be killed
#SBATCH --job-name DESeq2           # Job name
#SBATCH -o %j.out                       # File to which standard out will be written
#SBATCH -e %j.err               # File to which standard err will be written

module load gcc/6.2.0 R/3.5.1 hdf5/1.10.1


# This `for` loop will take the single-cell RNA-seq cluster ids as input and run the script for each of them on a different set of cores. The clusternames should be what is output on the seurat object from `seurat@ident` in the preparation for DESeq2 script.

for cluster_n in "clustername0" "clustername1" "clustername2" "clustername3" "clustername4" "clustername5" 

do

sbatch -p medium -t 3-12:00 -c 8 --mem 64G --job-name DEseq2 --wrap="Rscript sc_DESeq2_analysis_inner.R $cluster_n"

sleep 1 # wait 1 second between each job submission

done
