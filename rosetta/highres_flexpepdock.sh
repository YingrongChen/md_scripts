##!/bin/bash

#!/bin/csh

ROSETTA=/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin
# RELAX=/dors/meilerlab/home/chey120/chainA_chainA/Y39/relax
# HIGHRES=/dors/meilerlab/home/chey120/chainA_chainA/Y39/highres
DATA=/ssd1/chey120/workspace/chainA_chainA/chainA.db3
#rosetta command line
for file in `ls *.pdb`;
do
    name=${file##*/}
    name=${name%%_0001*} #remove _0001.pdb
    $ROSETTA/FlexPepDocking.bcl.linuxgccrelease -s $file -out:file:silent $name.silent -out:file:scorefile $name.sc -out:path:all $HIGHRES @highres_flag > $name.log &
done