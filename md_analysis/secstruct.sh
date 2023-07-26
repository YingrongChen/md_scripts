#bash secstruct.sh to_process

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
mask=":54-70"

# Loop over all proteins in the list
for dir in `cat ${LIST}`; do
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

if [ -f "${PROTEIN}_${Trial}_prod.offset_combine.nc" ]; then
#full distance
cat > secstruct.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.combine.nc
autoimage
secstruct ${mask} out ${PROTEIN}_${Trial}_secstruct.gnu sumout ${PROTEIN}_${Trial}_secstruct.dat
run
quit
EOF

cpptraj -i secstruct.in &

fi
done