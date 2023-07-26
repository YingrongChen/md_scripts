#!/bin/csh

ROSETTA=/ssd1/chey120/workspace/Rosetta/main/source/bin
RES=/ssd1/chey120/workspace/chainA_chainA/Y39_mut
OUT=/ssd1/chey120/workspace/chainA_chainA/Y39_mut
#rosetta command line
for resfile in $RES/*.resfile;
do
    for pdb in *.pdb;
    do
        name=${resfile##*/} 
        name=${name%.*}
        $ROSETTA/fixbb.bcl.linuxgccrelease -s $pdb -resfile $resfile -out:prefix ${name}_ -out:path:all $OUT
    done
done &