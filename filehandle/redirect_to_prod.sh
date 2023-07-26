#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/filehandle/redirect_to_prod.sh dir_list /dors/meilerlab/home/chey120/mhcii_asyn/MD/KEGVL_1h15

LIST=`readlink -e $1` # by readlink -e */ > dir_list
DIR=`readlink -e $2` #/dors/meilerlab/home/chey120/mhcii_asyn/MD/1bx2_prod

#rosetta command line
for i in `cat ${LIST}`; do
cd $i
PROTEIN=`basename $i`

#setup production run

for i in {0..4}; do
mkdir -p ${DIR}/${PROTEIN}/trial${i}
cp ${PROTEIN}_prod.0000.crd ${DIR}/${PROTEIN}/trial${i}
cp ${PROTEIN}.parm7 ${DIR}/${PROTEIN}/trial${i}
cp n_atoms.txt ${DIR}/${PROTEIN}/trial${i}
cp n_residues.txt ${DIR}/${PROTEIN}/trial${i}
cp stripped.${PROTEIN}.parm7 ${DIR}/${PROTEIN}/trial${i}

done
done &