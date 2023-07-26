# # "70@ND2190@O"
# awk '{ print $10 }' DNEAY_4x5w_hbond.dat > DNEAY_4x5w_hbond_190.dat
# awk '{ print $10 }' phos_DNEAY_4x5w_hbond.dat > phos_DNEAY_4x5w_hbond_190.dat
# awk '{ print $8 }' DNEAY_5ni9_hbond.dat > DNEAY_5ni9_hbond_190.dat
# awk '{ print $8 }' phos_DNEAY_5ni9_hbond.dat > phos_DNEAY_5ni9_hbond_190.dat

# # "113@OH190@N"
# awk '{ print $9 }' DNEAY_5ni9_hbond.dat > DNEAY_5ni9_hbond_190_113.dat
# awk '{ print $9 }' phos_DNEAY_5ni9_hbond.dat > phos_DNEAY_5ni9_hbond_190_113.dat

# # "144@NE1191@O"
# awk '{ print $11 }' DNEAY_4x5w_hbond.dat > DNEAY_4x5w_hbond_191.dat
# awk '{ print $11 }' phos_DNEAY_4x5w_hbond.dat > phos_DNEAY_4x5w_hbond_191.dat
# awk '{ print $10 }' DNEAY_5ni9_hbond.dat > DNEAY_5ni9_hbond_191.dat
# awk '{ print $10 }' phos_DNEAY_5ni9_hbond.dat > phos_DNEAY_5ni9_hbond_191.dat

# # "154@NH1188@O"    "154@NH2188@O"
# awk '{print $8<$9?$8:$9}' DNEAY_4x5w_hbond.dat > DNEAY_4x5w_hbond_188.dat 
# awk '{print $8<$9?$8:$9}' phos_DNEAY_4x5w_hbond.dat > phos_DNEAY_4x5w_hbond_188.dat 

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list

# Loop over all proteins in the list
for file in `cat ${LIST}`; do
file=`basename $file`
# awk '{ print $2 }' $file > temp_${file}
python /dors/meilerlab/home/chey120/chainA_chainA/scripts/md/dwell_time.py ${file}
done