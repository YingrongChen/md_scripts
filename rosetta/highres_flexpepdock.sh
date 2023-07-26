##!/bin/bash

#!/bin/csh

ROSETTA=rosetta/rosetta-3.13/main/source/bin
# RELAX=relax
# HIGHRES=highres
#rosetta command line
for file in `ls *.pdb`;
do
    name=${file##*/}
    name=${name%%_0001*} #remove _0001.pdb
    $ROSETTA/FlexPepDocking.bcl.linuxgccrelease -s $file -out:file:silent $name.silent -out:file:scorefile $name.sc -out:path:all $HIGHRES @highres_flag > $name.log &
done