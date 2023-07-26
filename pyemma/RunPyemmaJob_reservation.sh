#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=msm
#SBATCH --output=msm.%j.log
#SBATCH --reservation=md_simulation_brown_chen

#to submit:
#sbatch ../../scripts/pyemma/RunPyemmaJob_reservation.sh ../../scripts/pyemma/DNEAY_tica_build.py combine_DNEAY_5ni9.pdb combine_DNEAY_5ni9.nc
source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt
export PATH=/dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
echo $SLURM_JOB_NODELIST
echo $HOST $HOSTNAME
pv=`which python`
echo "Using $pv"

export HDF5_USE_FILE_LOCKING='FALSE'

# Input variables
WRAPPER=`readlink -e $1` # buildmsm.py
TOPFILE=`readlink -e $2` # pdb
TRAJFILE=`readlink -e $3` # traj
# PREFIX=$4 # output prefix3
# traj_size=$5

echo $WRAPPER $TOPFILE $TRAJFILE $PREFIX $traj_size
# if [ ! -f ${TRAJFILE/.nc/_image.nc} ]; then
# cat > image.in << EOF
# parm combine_DNEAY_5ni9.parm7
# trajin strip_DNEAY_5ni9.nc
# trajin strip_S129_DNEAY_5ni9.nc
# autoimage
# rms ToFirst out rmsd.dat @N,CA,C,O
# trajout combine_DNEAY_5ni9.nc
# run
# EOF
# cpptraj -i image.in 
# fi
# TRAJFILE=${TRAJFILE/.nc/_image.nc}

# python /dors/meilerlab/home/chey120/mhcii_asyn/scripts/pyemma/featurescore.py $TOPFILE $TRAJFILE $PREFIX
python $WRAPPER $TOPFILE $TRAJFILE $PREFIX $traj_size
