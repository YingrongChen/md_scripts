######################################################################
#Define and track “native” contacts as determined by a simple distance cut-off, i.e. any atoms which are closer than
#<cut> in the specified reference frame (the first frame if no reference specified) are considered a native contact. If
#one mask is provided, contacts are looked for within <mask1>; if two masks are provided, only contacts between
#atoms in <mask1> and atoms in <mask2> are looked for (useful for determining intermolecular contacts). By
#default only native contacts are tracked. This can be changed by specifying the savenonnative keyword. The time
#series for contacts can be saved using the series keyword; these can be further consolidated by residue using the
#resseries keyword. When using <resseries> the data set index is calculated as (r2 * nres) + r1 so that indices can
#be matched between native/non-native contact pairs. Non-native residue contact legends have an nn_ prefix.
#
#Above is extracted from Amber18  manual
#############################################################

#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md/nativecontacts.sh to_process
#source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
asyn=":180-196&!@H="
mhcii=":1-179&!@H="
cut="5.0"

# Loop over all proteins in the list
for dir in `cat ${LIST}`; do

cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

if [ -f "${PROTEIN}_${Trial}_prod.offset_combine.nc" ]; then

cat > nativecontacts.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.offset_combine.nc
autoimage
nativecontacts name NC1 \
    ${asyn} ${mhcii} \
    writecontacts ${PROTEIN}_${Trial}_nativecontacts.dat \
    distance ${cut} \
    byresidue out ${PROTEIN}_${Trial}_nativecontacts_residues.dat mindist maxdist \
    map mapout gnu \
    contactpdb ${PROTEIN}_${Trial}_nativecontacts_pdb.pdb \
    series seriesout ${PROTEIN}_${Trial}_nativecontacts_series.dat
run
lifetime NC1[NC] out ${PROTEIN}_${Trial}_nativecontacts_lifetime.dat
run
quit
EOF
cpptraj -i nativecontacts.in &

fi
done