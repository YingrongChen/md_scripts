#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md_analysis/mdout_nonTCR_process.sh to_process
source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

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

# bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/filehandle/rename.sh
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
asyn=":180-194"
mhcii=":1-179"
cut="5.0"
echo $Trial
echo $PROTEIN

######
# extract the data from the output files and plot
######
RESTART=`ls ${PROTEIN}_prod.*.crd | tail -n1 | awk -F. '{print $(NF-1)}'`
if [[ $RESTART == 0100 ]] && [ ! -f ${PROTEIN}_${Trial}_mhcii_rmsf.dat ]; then

######
# cpptraj: combine NC's to nc Binpos, calculate rmsf of peptide and hydrogen bond
######
cat > combine_nc.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_prod.0*.nc
autoimage
rms ToFirst out ${PROTEIN}_${Trial}_rmsd.dat @N,CA,C,O
atomicfluct ${asyn}@C,CA,N,O out ${PROTEIN}_${Trial}_asyn_rmsf.dat byres
atomicfluct ${mhcii}@C,CA,N,O out ${PROTEIN}_${Trial}_mhcii_rmsf.dat byres
hbond donormask ${asyn}@C,CA,N,O avgout ${PROTEIN}_${Trial}_backbond_hbt.dat nointramol
hbond acceptormask ${asyn}@C,CA,N,O avgout ${PROTEIN}_${Trial}_backbond_hbt.dat nointramol
nativecontacts name NC1 \
    ${asyn}&!@H= ${mhcii}&!@H= \
    writecontacts ${PROTEIN}_${Trial}_nativecontacts.dat \
    distance ${cut} \
    contactpdb ${PROTEIN}_${Trial}_nativecontacts_pdb.pdb
distance :62-72 :146-157 out ${PROTEIN}_${Trial}_mhciiwidth.dat
secstruct :62-72,146-157 out ${PROTEIN}_${Trial}_secstruct.gnu sumout ${PROTEIN}_${Trial}_secstruct.dat
trajout ${PROTEIN}_${Trial}_prod.offset_combine.nc offset 100
trajout ${PROTEIN}_${Trial}_prod.combine.nc
run
quit
EOF

#cpptraj -i combine_nc.in &

fi
done
)
# Wait for all background tasks to complete
#wait
(
for dir in `cat ${LIST} `; do
((n=n%${PROCESSORS})); ((n++==0)) && wait
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
if [ -f ${PROTEIN}_${Trial}_prod.offset_combine.nc ]; then
mkdir -p /dors/meilerlab/home/chey120/mhcii_asyn/local/${PROTEIN}
if [ ! -f "/dors/meilerlab/home/chey120/mhcii_asyn/local/${PROTEIN}/stripped.${PROTEIN}*.parm7" ]; then
cp stripped.${PROTEIN}*.parm7 /dors/meilerlab/home/chey120/mhcii_asyn/local/${PROTEIN}
fi
cp ${PROTEIN}_${Trial}_prod.offset_combine.nc /dors/meilerlab/home/chey120/mhcii_asyn/local/${PROTEIN} &
cp ${PROTEIN}_${Trial}_nativecontacts.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/nativecontacts_dat/mhcii_asyn
cp ${PROTEIN}_${Trial}_nativecontacts_pdb.pdb /dors/meilerlab/home/chey120/mhcii_asyn/Data/nativecontacts_pdb/mhcii_asyn
cp ${PROTEIN}_${Trial}_backbond_hbt.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/hbond_tcr
cp ${PROTEIN}_${Trial}_asyn_rmsf.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/rmsf/asyn
cp ${PROTEIN}_${Trial}_mhcii_rmsf.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/rmsf/mhcii
cp ${PROTEIN}_${Trial}_secstruct.gnu /dors/meilerlab/home/chey120/mhcii_asyn/Data/secstruct/
cp ${PROTEIN}_${Trial}_secstruct.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/secstruct/
cat ${PROTEIN}_${Trial}_mhciiwidth.dat >> /dors/meilerlab/home/chey120/mhcii_asyn/Data/special_distance/${PROTEIN}_width.dat
fi
done
) 


######
# cpptraj: Hbond with peptide
######
#specify backbone
# cat > hbond_native_contacts.in << EOF
# parm stripped.${PROTEIN}_cap_amber.parm7 
# trajin ${PROTEIN}_prod.combine.binpos
# hbond donormask :180-196 avgout avg-hbt.dat nointramol
# hbond acceptormask :180-196 avgout avg-hbt.dat nointramol
# EOF
# cpptraj -i hbond_native_contacts.in

######
# cpptraj: Hbond distribution
######
#specify backbone
# cat > hbond_native_contacts.in << EOF
# parm stripped.${PROTEIN}_cap_amber.parm7 
# trajin ${PROTEIN}_prod.combine.binpos
# atomicfluct :180-196@C,CA,N,O out peptide_rmsf.dat byres
# hbond donormask :180-196@C,CA,N,O avgout bbond_hbt.dat nointramol
# hbond acceptormask :180-196@C,CA,N,O avgout bbond_hbt.dat nointramol
# EOF
# cpptraj -i hbond_native_contacts.in

# ######
# # cpptraj: locate the structure with the lowest energy
# ######
# cat sum_EPTOT | awk '{if($2<min) {min=$2;print $1"   "min}}' > lowest_EPTOT
# lowest=`awk 'END {print $1}' lowest_EPTOT`
# lowest=${lowest%*.*} #drop decimal point
# frame=$((lowest/2500)) #ntwx=2500r
# file=`grep $lowest *.out | awk '{print $1}'`
# file=${file%*.out:} #drop .out
# cat > extract_${frame}.in << EOF
# parm stripped.${PROTEIN}_cap_amber.parm7 
# trajin ${file}.nc $frame $frame
# trajout lowest_energy_struct.pdb pdb
# EOF
# cpptraj -i extract_${frame}.in

# ######
# # cpptraj: rmsd to the lowest energy structure
# ######
# cat > rmsd_to_${frame}.in << EOF
# parm stripped.${PROTEIN}_cap_amber.parm7 
# trajin ${PROTEIN}_prod.0*.nc
# reference lowest_energy_struct.pdb
# rms reference out rmsd_to_${frame}.dat @N,CA,C time 1.0 
# EOF
# cpptraj -i rmsd_to_${frame}.in

    # for run in `seq 0 1 4`; do
    #     mkdir trial${run}
    #     cd trial${run}
    #     cp ../n_atoms.txt .
    #     cp ../n_residues.txt .
    #     cp ../${PROTEIN}.parm7 .
    #     cp ../${PROTEIN}_prod.0000.crd .
    #     cp ../stripped.${PROTEIN}.parm7 .
    #     rm -r trial*/
    #     cd ../
    # done
