#bash crd_out.sh to_process

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
PROCESSORS=10			# Number of parallel processes

(
for dir in `cat ${LIST} `; do
((n=n%${PROCESSORS})); ((n++==0)) && wait
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
echo $Trial $PROTEIN
RESTART=`ls ${PROTEIN}_prod.*.crd | tail -n1 | awk -F. '{print $(NF-1)}'`
if [[ $RESTART == 0100 ]] && [ -f ${PROTEIN}_${Trial}_chainA_rmsf.dat ]; then
cat > crd.in << EOF
parm ${PROTEIN}*.parm7 
trajin ${PROTEIN}_prod.0*.crd
autoimage
trajout ${PROTEIN}_${Trial}_crd.nc
run
EOF

cpptraj -i crd.in &

cat > radial.in << EOF
parm ${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_crd.nc
radial out 200mM-Na.dat 0.1 15.0 :189 :Na+
radial out 200mM-Cl.dat 0.1 15.0 :189 :Cl-
run
EOF

cpptraj -i radial.in &
fi
done
)