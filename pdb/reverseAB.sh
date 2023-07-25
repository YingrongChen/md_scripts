for PDB in *5ni9*.pdb;
do
    tag=`basename $PDB .pdb`
    # Run
    awk '$5=="A"{print}' $PDB> ${tag}_r.pdb
    awk '$5=="B"{print}' $PDB>> ${tag}_r.pdb
    awk '$5=="C"{print}' $PDB>> ${tag}_r.pdb
    tail -3 $PDB>> ${tag}_r.pdb
done