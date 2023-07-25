#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=secstruct
#SBATCH --output=secstruct.%j.log
#SBATCH --reservation=md_simulation_brown_chen

#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md_analysis/secstruct.sh to_process
source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

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
atominfo :189
#N-CA-CB-CG	
dihedral ${PROTEIN}_${Trial}_tyrchi1 :189@N :189@CA :189@CB :189@OG out ${PROTEIN}_${Trial}_tyrchi1.dat
#CA-CB-CG-CD1
#dihedral ${PROTEIN}_${Trial}_tyrchi2 :189@CA :189@CB :189@CG :189@CD1 out ${PROTEIN}_${Trial}_tyrchi2.dat
run
quit
EOF
cpptraj -i secstruct_5ni9.in &
fi
done