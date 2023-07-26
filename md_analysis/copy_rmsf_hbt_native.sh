# bash /copy_rmsf_hbt_native.sh to_process
LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
for dir in `cat ${LIST}`; do
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
echo $PROTEIN
# # if [ $TCR == "TRUE" ]; then
# cat ${PROTEIN}_${Trial}_chainAwidth.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_width.dat
# cat ${PROTEIN}_${Trial}_mhcdTCR.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_mhcdTCR.dat
# cat ${PROTEIN}_${Trial}_chainAdTCR.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_chainAdTCR.dat
# cat ${PROTEIN}_${Trial}_94d279.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_94d279.dat
# cat ${PROTEIN}_${Trial}_Bd286.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_Bd286.dat
# cat ${PROTEIN}_${Trial}_tcrwidth1.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_tcrwidth1.dat
# cat ${PROTEIN}_${Trial}_tcrwidth2.dat >> /dors/meilerlab/home/chey120/chainA_chainA/Data/special_distance/${PROTEIN}_tcrwidth2.dat
cp ${PROTEIN}_${Trial}_tyrchi1.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/dihedral
# cp ${PROTEIN}_${Trial}_tyrchi2.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/dihedral

# cp ${PROTEIN}_${Trial}_409hb203.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/hbond_distance
# cp ${PROTEIN}_${Trial}_380hb64.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/hbond_distance

cp ${PROTEIN}_${Trial}_secstruct.gnu /dors/meilerlab/home/chey120/chainA_chainA/Data/secstruct/
cp ${PROTEIN}_${Trial}_secstruct.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/secstruct/

# cp ${PROTEIN}_${Trial}_prod.offset_combine.nc /dors/meilerlab/home/chey120/chainA_chainA/local/${PROTEIN} &
# cp ${PROTEIN}_${Trial}_chainA_tcr_nativecontacts.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/nativecontacts_dat/tcr_chainA
# cp ${PROTEIN}_${Trial}_chainA_tcr_nativecontacts.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/nativecontacts_dat/tcr_chainA
# cp ${PROTEIN}_${Trial}_backbond_hbt.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/hbond_tcr
# cp ${PROTEIN}_${Trial}_chainA_rmsf.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/rmsf/chainA
# cp ${PROTEIN}_${Trial}_chainA_rmsf.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/rmsf/chainA
# cp ${PROTEIN}_${Trial}_tcr_rmsf.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/rmsf/tcr
# cp ${PROTEIN}_${Trial}_nativecontacts_pdb.pdb /dors/meilerlab/home/chey120/chainA_chainA/Data/nativecontacts_pdb/chainA_chainA
# cp ${PROTEIN}_${Trial}_prod.offset_combine.nc /dors/meilerlab/home/chey120/chainA_chainA/local/${PROTEIN} &
# cp ${PROTEIN}_${Trial}_nativecontacts.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/nativecontacts_dat/chainA_chainA
# cp ${PROTEIN}_${Trial}_backbond_hbt.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/hbond_tcr
# cp ${PROTEIN}_${Trial}_peptide_rmsf.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/rmsf/chainA
# cp ${PROTEIN}_${Trial}_chainA_rmsf.dat /dors/meilerlab/home/chey120/chainA_chainA/Data/rmsf/chainA
# fi

    # sed "s/$/\t$PROTEIN/" peptide_rmsf.dat > named_peptide_rmsf.dat
    # cat named_peptide_rmsf.dat >> ../../${Trial}_peptide_rmsf.dat

done