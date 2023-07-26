#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=dihedral
#SBATCH --output=dihedral.%j.log
#SBATCH --reservation=md_simulation_brown_chen

#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md_analysis/dihedral.sh to_process
source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list

# Loop over all proteins in the list
for dir in `cat ${LIST}`; do

cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

if [ -f "*KEGVL*prod.offset_combine.nc" ]; then
#full distance
cat > dihedral.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.combine.nc
autoimage
atominfo :184
#N-CA-CB-CG	
dihedral ${PROTEIN}_${Trial}_tyrchi1 :184@N :184@CA :184@CB :184@CG out ${PROTEIN}_${Trial}_tyrchi1.dat
#CA-CB-CG-CD1
dihedral ${PROTEIN}_${Trial}_tyrchi2 :184@CA :184@CB :184@CG :184@CD1 out ${PROTEIN}_${Trial}_tyrchi2.dat
run
quit
EOF
cpptraj -i dihedral.in &
fi
if [ -f "*DNEAY_5ni9*prod.offset_combine.nc" ]; then
#full distance
cat > dihedral.in << EOF
parm stripped.${PROTEIN}*.parm7 
trajin ${PROTEIN}_${Trial}_prod.combine.nc
autoimage
atominfo :184
#N-CA-CB-CG	
dihedral ${PROTEIN}_${Trial}_tyrchi1 :184@N :184@CA :184@CB :184@CG out ${PROTEIN}_${Trial}_tyrchi1.dat
#CA-CB-CG-CD1
dihedral ${PROTEIN}_${Trial}_tyrchi2 :184@CA :184@CB :184@CG :184@CD1 out ${PROTEIN}_${Trial}_tyrchi2.dat
run
quit
EOF
cpptraj -i dihedral.in &
fi
done