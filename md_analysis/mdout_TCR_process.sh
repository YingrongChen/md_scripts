#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=mdout
#SBATCH --output=mdout.%j.log
#SBATCH --mail-user=ych2429@emory.edu
#SBATCH --mail-type=END,FAIL

source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt
# bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md_analysis/mdout_TCR_process.sh to_process

# export PATH=/dors/meilerlab/home/chey120/mhcii_asyn/scripts/md:${PATH}
# chmod u+x /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md/mdout_extract.perl
# module load GCC OpenMPI R

#find . -name prod_47381594.log

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list

PROCESSORS=10			# Number of parallel processes

# Parallel process for-loop without using xargs
(
for dir in `cat ${LIST} `; do
((n=n%${PROCESSORS})); ((n++==0)) && wait
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
asynend=$(cat n_residues.txt)
asynstart=$(($asynend-14))
mhciiend=$(($asynstart-1))
mhciistart=$(($mhciiend-177))
tcrend=$(($mhciistart-1))
asyn=":${asynstart}-${asynend}"
asyn=":${mhciistart}-${mhciiend}"
asyn=":1-${tcrend}"
cut="5.0"
ptr=$(($asynend-9))
# cp /dors/meilerlab/home/chey120/mhcii_asyn/MD/init/TCR/$PROTEIN/stripped.${PROTEIN}.parm7 .
# ######
# # extract the data from the output files and plot
# ######
RESTART=`ls ${PROTEIN}_prod.*.crd | tail -n1 | awk -F. '{print $(NF-1)}'`
# if [ ! -f ${PROTEIN}_${Trial}_mhcii_rmsf.dat ] && [[ $RESTART == 0100 ]] ; then

# ######
# # cpptraj: combine NC's to nc Binpos, calculate rmsf of peptide and hydrogen bond
# ######
cat > combine_nc.in << EOF
parm stripped.${PROTEIN}*.parm7
trajin ${PROTEIN}_prod.0*.nc
autoimage
rms ToFirst out ${PROTEIN}_${Trial}_rmsd.dat @N,CA,C,O
trajout ${PROTEIN}_${Trial}_prod.offset_combine.nc offset 100
trajout ${PROTEIN}_${Trial}_prod.combine.nc
atomicfluct ${asyn}@C,CA,N,O out ${PROTEIN}_${Trial}_asyn_rmsf.dat byres
atomicfluct ${mhcii}@C,CA,N,O out ${PROTEIN}_${Trial}_mhcii_rmsf.dat byres
atomicfluct ${tcr}@C,CA,N,O out ${PROTEIN}_${Trial}_tcr_rmsf.dat byres
hbond donormask ${asyn}@C,CA,N,O avgout ${PROTEIN}_${Trial}_backbond_hbt.dat nointramol
hbond acceptormask ${asyn}@C,CA,N,O avgout ${PROTEIN}_${Trial}_backbond_hbt.dat nointramol
nativecontacts name NC1 \
    ${tcr}&!@H= ${mhcii}&!@H= \
    writecontacts ${PROTEIN}_${Trial}_mhcii_tcr_nativecontacts.dat \
    distance ${cut} \
    contactpdb ${PROTEIN}_${Trial}_mhcii_tcr_nativecontacts.pdb    
nativecontacts name NC2 \
    ${tcr}&!@H= ${asyn}&!@H= \
    writecontacts ${PROTEIN}_${Trial}_asyn_tcr_nativecontacts.dat \
    distance ${cut} \
    contactpdb ${PROTEIN}_${Trial}_asyn_tcr_nativecontacts.pdb 
#N-CA-CB-CG	
dihedral ${PROTEIN}_${Trial}_tyrchi1 :${ptr}@N :${ptr}@CA :${ptr}@CB :${ptr}@CG out ${PROTEIN}_${Trial}_tyrchi1.dat
#CA-CB-CG-CD1
dihedral ${PROTEIN}_${Trial}_tyrchi2 :${ptr}@CA :${ptr}@CB :${ptr}@CG :${ptr}@CD1 out ${PROTEIN}_${Trial}_tyrchi2.dat
run
quit
EOF
cpptraj -i combine_nc.in &
echo $dir
# fi
done
)

# for dir in `cat ${LIST}`; do
# cd $dir
# Trial=`basename $dir`
# PROTEIN=${dir%*/trial*} #remove /trial0
# PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
# if [ -f ${PROTEIN}_${Trial}_prod.offset_combine.nc ]; then
# echo $PROTEIN $Trial
# cp ${PROTEIN}_${Trial}_prod.offset_combine.nc /dors/meilerlab/home/chey120/mhcii_asyn/local/${PROTEIN} &
# cp ${PROTEIN}_${Trial}_mhcii_tcr_nativecontacts.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/nativecontacts_dat/tcr_mhcii
# cp ${PROTEIN}_${Trial}_asyn_tcr_nativecontacts.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/nativecontacts_dat/tcr_asyn
# cp ${PROTEIN}_${Trial}_backbond_hbt.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/hbond_tcr
# cp ${PROTEIN}_${Trial}_asyn_rmsf.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/rmsf/asyn
# cp ${PROTEIN}_${Trial}_mhcii_rmsf.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/rmsf/mhcii
# cp ${PROTEIN}_${Trial}_tcr_rmsf.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/rmsf/tcr
# fi
# done