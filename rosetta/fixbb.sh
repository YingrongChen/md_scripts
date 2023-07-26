#!/bin/csh

ROSETTA=Rosetta/main/source/bin
RES=
OUT=
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