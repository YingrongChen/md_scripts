#!/bin/bash

#to submit:
#sbatch ../../scripts/pyemma/RunPyemmaJob.sh backbone_nontcr.pdb backbone_nontcr_prod.nc new_
export PATH=miniconda3/envs/pyemma/bin:$PATH
echo $SLURM_JOB_NODELIST
echo $HOST $HOSTNAME
pv=`which python`
echo "Using $pv"

export HDF5_USE_FILE_LOCKING='FALSE'

PROTEIN=$(basename "$(pwd)")
if [[ -f bb_${PROTEIN}*.combine.nc ]]; then
cp MD/prod/$PROTEIN*/trial?/bb_${PROTEIN}_trial?_4x5w.combine.nc
fi
# Input variables
WRAPPER=`readlink -e $1` # buildmsm.py
TOPFILE=`readlink -e $2` # pdb
TRAJFILE=`readlink -e $3` # traj
PREFIX=$4 # output prefix3
traj_size=$5

echo $WRAPPER $TOPFILE $TRAJFILE $PREFIX $traj_size
# if [ ! -f ${TRAJFILE/.nc/_image.nc} ]; then
# cat > image.in << EOF
# parm bb.bb.stripped.S129_DNEAY_4x5w.parm7
# trajin $TRAJFILE
# autoimage
# rms ToFirst out rmsd.dat @N,CA,C,O
# trajout ${TRAJFILE/.nc/_image.nc}
# run
# EOF
# cpptraj -i image.in 
# fi
# TRAJFILE=${TRAJFILE/.nc/_image.nc}

# python pyemma/featurescore.py $TOPFILE $TRAJFILE $PREFIX
python $WRAPPER $TOPFILE $TRAJFILE $PREFIX $traj_size
