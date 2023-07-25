#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=cpptraj
#SBATCH --output=cpptraj.%j.log
#SBATCH --reservation=md_simulation_brown_chen

#to submit:
# sbatch ../scripts/pyemma/cpptraj.sh

source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

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