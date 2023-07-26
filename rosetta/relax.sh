#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=300G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=relax
#SBATCH --output=relax.%j.log
#SBATCH --reservation=md_simulation_brown_chen

#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/rosetta/relax.sh

# for file in *-1.pdb; do
#     cat DNEAY_5ni9_1.pdb >> "$file"
# done

ROSETTA=/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/dors/meilerlab/apps/Linux2/x86_64/lib64/:/dors/meilerlab/apps/Linux2/x86_64/lib/
# DIR_LIST=$1
#rosetta command line
# for DIR in `cat ${DIR_LIST}`;
# do
# cd $DIR
for PDB in `ls *.pdb`;
do
name=${PDB##*/} 
name=${PDB%.*}

$ROSETTA/relax.linuxgccrelease \
-out:prefix ${name}_ \
-in:file:s $PDB \
-relax:constrain_relax_to_start_coords true \
-coord_constrain_sidechains false \
-ramp_constraints false \
-score:weights ref2015 > ${name}.log &

done