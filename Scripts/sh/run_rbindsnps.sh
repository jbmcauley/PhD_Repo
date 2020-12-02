#!/bin/bash

#$ -V
#$ -cwd
#$ -N rbind_snpobjects
#$ -o o_files/
#$ -e e_files/

#activate genedrop conda environment
source /ceph/software/conda/etc/profile.d/conda.sh &&
conda activate jmcauley_rbindsnpobjs &&

# make directory in /scratch for this job
SCRATCH=/scratch/$USER/rbindsnpobjs/ &&
mkdir -p $SCRATCH &&

# run command
Rscript scripts/example_gene_dropping_script_dummy.R $SCRATCH &&

# move the result file back to the directory the job was submitted from using rsync
rsync -av $SCRATCH /data/johnston/sparrows/analyses/results/ &&

# remove the directory ceated for this job to clear up any unwanted files
#rm -rf $SCRATCH 

#deactivate conda environment
conda deactivate