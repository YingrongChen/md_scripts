#!/bin/bash
#SBATCH --account=csb_gpu_acc
#SBATCH --job-name=mmpbsa
#SBATCH --output=mmpbsa_%j.log
#SBATCH --partition=turing
#SBATCH --gres-flags=enforce-binding
#SBATCH --gres=gpu:1
#SBATCH --mem=6G
#SBATCH --time=5-00:00:00
#SBATCH --mail-user=ych2429@emory.edu
#SBATCH --mail-type=END,FAIL

source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

# LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
# for traj in `cat ${LIST}`; do
# Trial=`basename $dir`
# PROTEIN=${dir%*/trial*} #remove /trial0
# ante-MMPBSA.py -p stripped.TCR_DNEAY.parm7 -c com.TCR_DNEAY.parm7 -r tcr.TCR_DNEAY.parm7 -l pMHC.TCR_DNEAY.parm7 -n :222-415

cat > mmpbsa.in << EOF
Input file for running PB and GB
&general
   endframe=50, verbose=1,
#  entropy=1,
/
&gb
  igb=2, saltcon=0.100
/
&pb
  istrng=0.100,
/
EOF

$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa.in \
    -o TCR_DNEAY_MMPBSA.dat \
    -cp stripped.TCR_DNEAY.parm7 \
    -rp tcr.TCR_DNEAY.parm7 \
    -lp pMHC.TCR_DNEAY.parm7 \
    -y TCR_DNEAY_prod.0100.nc

