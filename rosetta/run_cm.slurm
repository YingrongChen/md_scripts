#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=rosettacm
#SBATCH --output=rosettacm.log
#SBATCH --mail-user=ych2429@emory.edu
#SBATCH --mail-type=END,FAIL

seed=$SLURM_JOBID

/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.linuxgccrelease @rosetta_cm.options