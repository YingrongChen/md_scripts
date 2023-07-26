#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md_analysis/secstruct.sh to_process

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
#full distance
cat > secstruct_5ni9.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.combine.nc
autoimage
#TCR-pMHC
secstruct :54-70 out ${PROTEIN}_${Trial}_secstruct.gnu sumout ${PROTEIN}_${Trial}_secstruct.dat
#secstruct :275-291 out ${PROTEIN}_${Trial}_secstruct.gnu sumout ${PROTEIN}_${Trial}_secstruct.dat
run
quit
EOF
cpptraj -i secstruct_5ni9.in &
fi
done