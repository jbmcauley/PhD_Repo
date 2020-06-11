#!/bin/sh
#$ -N test              
#$ -cwd
#$ -l h_rt=00:30:00
#$ -V
#$ -l h_vmem=5200M

. /etc/profile.d/modules.sh
module load R

R CMD BATCH test.R
