#!/bin/bash

ROSETTA=rosetta/rosetta-3.13/main/source/bin/

readlink -e *.pdb > pdb.list
file=`awk NR==${SLURM_ARRAY_TASK_ID} pdb.list`
name=${file##*/}
name=${name%%_0001*} #remove _0001.pdb
$ROSETTA/FlexPepDocking.default.linuxgccrelease -s $file -out:file:silent $name.silent -out:file:scorefile $name.sc @highres_flag > $name.log
