#!/bin/bash
#SBATCH --account=p_meiler
#SBATCH --job-name=highres_5ni9
#SBATCH --output=highres_%j.log
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --array=1-64

ROSETTA=/gpfs51/dors2/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/

readlink -e *.pdb > pdb.list
file=`awk NR==${SLURM_ARRAY_TASK_ID} pdb.list`
name=${file##*/}
name=${name%%_0001*} #remove _0001.pdb
$ROSETTA/FlexPepDocking.default.linuxgccrelease -s $file -out:file:silent $name.silent -out:file:scorefile $name.sc @highres_flag > $name.log
