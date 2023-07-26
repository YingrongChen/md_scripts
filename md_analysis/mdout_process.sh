######
# cpptraj: combine NC's to a combine and a offseted versions, calculate RMSF, hydrogen bond, native contacts between two chains
# bash md_analysis/mdout_nonTCR_process.sh to_process
######
LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
PROCESSORS=10			# Number of parallel processes
chainA=":180-194"
chainB=":1-179"
(
for dir in `cat ${LIST} `; do
((n=n%${PROCESSORS})); ((n++==0)) && wait
cd $dir

Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
echo $Trial
echo $PROTEIN

RESTART=`ls ${PROTEIN}_prod.*.crd | tail -n1 | awk -F. '{print $(NF-1)}'`
if [[ $RESTART == 0100 ]] && [ ! -f ${PROTEIN}_${Trial}_chainB_rmsf.dat ]; then
cat > combine_nc.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_prod.0*.nc
autoimage
rms ToFirst out ${PROTEIN}_${Trial}_rmsd.dat @N,CA,C,O
atomicfluct ${chainA}@C,CA,N,O out ${PROTEIN}_${Trial}_chainA_rmsf.dat byres
atomicfluct ${chainB}@C,CA,N,O out ${PROTEIN}_${Trial}_chainB_rmsf.dat byres
hbond donormask ${chainA}@C,CA,N,O avgout ${PROTEIN}_${Trial}_backbond_hbt.dat nointramol
hbond acceptormask ${chainA}@C,CA,N,O avgout ${PROTEIN}_${Trial}_backbond_hbt.dat nointramol
nativecontacts name NC1 \
    ${chainA}&!@H= ${chainB}&!@H= \
    writecontacts ${PROTEIN}_${Trial}_nativecontacts.dat \
    contactpdb ${PROTEIN}_${Trial}_nativecontacts_pdb.pdb
trajout ${PROTEIN}_${Trial}_prod.offset_combine.nc offset 100
trajout ${PROTEIN}_${Trial}_prod.combine.nc
run
quit
EOF
cpptraj -i combine_nc.in &
fi
done
)
# Wait for all background tasks to complete
wait
(
for dir in `cat ${LIST} `; do
((n=n%${PROCESSORS})); ((n++==0)) && wait
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
if [ -f ${PROTEIN}_${Trial}_prod.offset_combine.nc ]; then
mkdir -p /local/${PROTEIN}
cp stripped.${PROTEIN}*.parm7 /local/${PROTEIN}
cp ${PROTEIN}_${Trial}_prod.offset_combine.nc /local/${PROTEIN} &
cp ${PROTEIN}_${Trial}_nativecontacts.dat /Data/nativecontacts_dat/chainB_chainA
cp ${PROTEIN}_${Trial}_nativecontacts_pdb.pdb /Data/nativecontacts_pdb/chainB_chainA
cp ${PROTEIN}_${Trial}_backbond_hbt.dat /Data/hbond_tcr
cp ${PROTEIN}_${Trial}_chainA_rmsf.dat /Data/rmsf/chainA
cp ${PROTEIN}_${Trial}_chainB_rmsf.dat /Data/rmsf/chainB
fi
done
) 