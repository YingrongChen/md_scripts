#!/bin/bash
#to submit:
#sbatch pyemma/RunPyemmaJob_reservation.sh pyemma/projectmsm.py [pdb] [nc]
export PATH=miniconda3/envs/pyemma/bin:$PATH
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

# python pyemma/featurescore.py $TOPFILE $TRAJFILE $PREFIX
python $WRAPPER $TOPFILE $TRAJFILE $PREFIX $traj_size
