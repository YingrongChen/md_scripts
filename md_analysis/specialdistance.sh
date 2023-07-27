#!/bin/bash
LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
chainA=":180-196&!@H="
chainA=":1-179&!@H="
cut="5.0"

# Loop over all proteins in the list
for dir in `cat ${LIST}`; do

cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

if [ -f "${PROTEIN}_${Trial}_prod.offset_combine.nc" ]; then

# cat > hbonddistance_4x5w.in << EOF
# parm stripped.${PROTEIN}*.parm7 
# trajin ${PROTEIN}_${Trial}_prod.combine.nc
# autoimage

# distance :409@O3P :203@NH2 out ${PROTEIN}_${Trial}_409hb203.dat
# distance :409@O2P :203@NE out ${PROTEIN}_${Trial}_409hb203.dat
# distance :409@O2P :203@NH2 out ${PROTEIN}_${Trial}_409hb203.dat
# distance :409@O1P :203@NH2 out ${PROTEIN}_${Trial}_409hb203.dat
# distance :409@O1P :203@NE out ${PROTEIN}_${Trial}_409hb203.dat

# distance :380@OD1 :64@NH1 out ${PROTEIN}_${Trial}_380hb64.dat
# distance :380@OD1 :64@NH2 out ${PROTEIN}_${Trial}_380hb64.dat
# distance :380@OD2 :64@NH1 out ${PROTEIN}_${Trial}_380hb64.dat
# distance :380@OD2 :64@NH2 out ${PROTEIN}_${Trial}_380hb64.dat

# run
# quit
# EOF
# cpptraj -i hbonddistance_4x5w.in &

# cat > hbonddistance_5ni9.in << EOF
# parm stripped.${PROTEIN}*.parm7 
# trajin ${PROTEIN}_${Trial}_prod.combine.nc
# autoimage
# #TCR-pMHC
# distance :188@O3P :397@NH2 out ${PROTEIN}_${Trial}_188hb397.dat
# distance :188@O2P :397@NE out ${PROTEIN}_${Trial}_188hb397.dat
# distance :188@O2P :397@NH2 out ${PROTEIN}_${Trial}_188hb397.dat
# distance :188@O1P :397@NH2 out ${PROTEIN}_${Trial}_188hb397.dat
# distance :188@O1P :397@NE out ${PROTEIN}_${Trial}_188hb397.dat

# distance :159@OD1 :258@NH1 out ${PROTEIN}_${Trial}_159hb258.dat
# distance :159@OD1 :258@NH2 out ${PROTEIN}_${Trial}_159hb258.dat
# distance :159@OD2 :258@NH1 out ${PROTEIN}_${Trial}_159hb258.dat
# distance :159@OD2 :258@NH2 out ${PROTEIN}_${Trial}_159hb258.dat

# run
# quit
# EOF
# cpptraj -i hbonddistance_5ni9.in

# cat > specialdistance_4x5w.in << EOF
# parm stripped.${PROTEIN}*.parm7 
# trajin ${PROTEIN}_${Trial}_prod.combine.nc
# autoimage
# #TCR-pMHC
# distance :48@CB,N,CA,O,C :377@CB out ${PROTEIN}_${Trial}_48d377_full.dat
# distance :137 :286-287 out ${PROTEIN}_${Trial}_137d286.dat
# distance :47-49 :157-159 out ${PROTEIN}_${Trial}_tcrwidth2.dat
# distance :47-49 :203-208 out ${PROTEIN}_${Trial}_tcrwidth1.dat
# distance :64@CZ :380@CG out ${PROTEIN}_${Trial}_64d380.dat
# distance :283-293 :367-378 out ${PROTEIN}_${Trial}_chainAwidth.dat
# distance :401-415 :1-221 out ${PROTEIN}_${Trial}_chainAdTCR.dat
# distance :222-400 :1-221 out ${PROTEIN}_${Trial}_mhcdTCR.dat
# distance :94 :279 out ${PROTEIN}_${Trial}_94d279.dat
# distance :111-221 :286-287 out ${PROTEIN}_${Trial}_Bd286.dat
# run
# quit
# EOF
# cpptraj -i specialdistance_4x5w.in

cat > specialdistance_5ni9.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.combine.nc
autoimage
#TCR-pMHC
distance :241-243 :377-379 out ${PROTEIN}_${Trial}_tcrwidth2.dat
distance :241-243 :203-208 out ${PROTEIN}_${Trial}_tcrwidth1.dat
distance :62-72 :146-157 out ${PROTEIN}_${Trial}_chainAwidth.dat
distance :180-194 :195-415 out ${PROTEIN}_${Trial}_chainAdTCR.dat
distance :1-179 :195-415 out ${PROTEIN}_${Trial}_mhcdTCR.dat
distance :288 :58 out ${PROTEIN}_${Trial}_94d279.dat
distance :305-415 :65-66 out ${PROTEIN}_${Trial}_Bd286.dat
run
quit
EOF
cpptraj -i specialdistance_5ni9.in

# cat > hbonddistance_1bx2.in << EOF
# parm stripped.${PROTEIN}*.parm7 
# trajin ${PROTEIN}_${Trial}_prod.combine.nc
# autoimage
# #TCR-pMHC
# distance :187@O3P :396@NH2 out ${PROTEIN}_${Trial}_187hb396.dat
# distance :187@O2P :396@NE out ${PROTEIN}_${Trial}_187hb396.dat
# distance :187@O2P :396@NH2 out ${PROTEIN}_${Trial}_187hb396.dat
# distance :187@O1P :396@NH2 out ${PROTEIN}_${Trial}_187hb396.dat
# distance :187@O1P :396@NE out ${PROTEIN}_${Trial}_187hb396.dat

# distance :158@OD1 :257@NH1 out ${PROTEIN}_${Trial}_158hb257.dat
# distance :158@OD1 :257@NH2 out ${PROTEIN}_${Trial}_158hb257.dat
# distance :158@OD2 :257@NH1 out ${PROTEIN}_${Trial}_158hb257.dat
# distance :158@OD2 :257@NH2 out ${PROTEIN}_${Trial}_158hb257.dat

# run
# quit
# EOF
# cpptraj -i hbonddistance_1bx2.in

fi
done