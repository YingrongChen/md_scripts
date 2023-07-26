#!/bin/bash

#bash rosetta/relax.sh

# for file in *-1.pdb; do
#     cat DNEAY_5ni9_1.pdb >> "$file"
# done

ROSETTA=rosetta/rosetta-3.13/main/source/bin

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