#!/bin/bash

#to submit:
# sbatch ../scripts/pyemma/cpptraj.sh


# protein_list=("DNEAY_5ni9" "S129_DNEAY_5ni9" "Y125_DNEAY_5ni9")

# for protein in "${protein_list[@]}"
# do
# cat > image.in << EOF
# parm stripped.${protein}*.parm7
# trajin ${protein}_trial?_prod.combine.nc
# strip @H*,P,O1P,O2P,O3P|:185@OH
# trajout strip_${protein}.nc
# trajout combine_DNEAY_5ni9.nc append
# run
# EOF
# cpptraj -i image.in 
# done