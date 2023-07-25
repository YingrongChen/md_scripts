#!/bin/bash

# visit each dir with multiple structures file and one score file
# compute rmsd with respect to the lowest energy geometry

ROSETTA=/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin
ANALYSIS=/dors/meilerlab/home/chey120/mhcii_asyn/Y39/analysis
SCRIPTS=/dors/meilerlab/home/chey120/mhcii_asyn/scripts


for dir in $(ls -d */); do

    name=${dir%%/*} 
    cd $name
    echo $name
    if [ ! -f rmsd_* ]; then
        lowest=$(awk 'FNR == 2{print $(NF)}' sort_$name.sc)

        $ROSETTA/extract_pdbs.linuxgccrelease \
            -in::file::silent $name.silent \
            -in::file::silent_struct_type binary \
            -in:file:tags $lowest

        $ROSETTA/rosetta_scripts.default.linuxgccrelease \
            -in::file::silent $name.silent \
            -in::file::silent_struct_type binary \
            -in:file:native $lowest.pdb \
            -parser:protocol  $SCRIPTS/ddg_rmsd.xml \
            -out:file:scorefile rmsd_$name.sc \
            -mute protocols.rosetta_scripts.ParsedProtocol.REPORT > $name.log

        sed -i '/^SEQUENCE/d' rmsd_$name.sc
        #Rscript $SCRIPTS/score_vs_rmsd.R rmsd_$name.sc
    else
        echo "rmsd computed"
    fi
    cd $ANALYSIS
done &

#/dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease -in:file:native 4gg6.pdb -parser:protocol /dors/meilerlab/home/chey120/mhcii_asyn/scripts/rosetta/rmsd_lowest_ddg_hbond.xml -out:file:scorefile rmsd.sc
