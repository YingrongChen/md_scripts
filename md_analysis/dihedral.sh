#!/bin/bash

#bash /dihedral.sh to_process

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list

# Loop over all proteins in the list
for dir in `cat ${LIST}`; do
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
cat > dihedral.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.combine.nc
autoimage
atominfo :184
#N-CA-CB-CG	
dihedral ${PROTEIN}_${Trial}_chi1 :184@N :184@CA :184@CB :184@CG out ${PROTEIN}_${Trial}_chi1.dat
#CA-CB-CG-CD1
dihedral ${PROTEIN}_${Trial}_chi2 :184@CA :184@CB :184@CG :184@CD1 out ${PROTEIN}_${Trial}_chi2.dat
run
quit
EOF
cpptraj -i dihedral.in &
done