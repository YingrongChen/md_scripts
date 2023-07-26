#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=cpptraj
#SBATCH --output=cpptraj.%j.log
#SBATCH --mail-user=ych2429@emory.edu
#SBATCH --mail-type=END,FAIL

source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
END_DIR=/dors/meilerlab/home/chey120/mhcii_asyn/commontraj
for dir in `cat ${LIST}`; do
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

if [[ !($PROTEIN =~ .*TCR.*) && $PROTEIN =~ .*4x5w.* ]] ; then
cat > stripparm.in << EOF
parm stripped.${PROTEIN}*.parm7
trajin ${PROTEIN}_${Trial}_prod.combine.nc
strip !(@C,CA,N,O)|:85,180,196 outprefix bb 
trajout bb_${PROTEIN}_${Trial}_4x5w.combine.nc
trajout ${END_DIR}/bb_nontcr_4x5w.nc append
trajout ${END_DIR}/bb_nontcr.nc append
EOF
cpptraj -i stripparm.in > stripparm.log

elif [[ !($PROTEIN =~ .*TCR.*) && $PROTEIN =~ .*5ni9.* ]] ; then
cat > stripparm.in << EOF
parm stripped.${PROTEIN}*.parm7
trajin ${PROTEIN}_${Trial}_prod.combine.nc
strip !(@C,CA,N,O)|:85,180,196 outprefix bb 
trajout backbone_${PROTEIN}_${Trial}_5ni9.combine.nc
trajout ${END_DIR}/bb_nontcr_5ni9.nc append
trajout ${END_DIR}/bb_nontcr.nc append
EOF
cpptraj -i stripparm.in > stripparm.log

elif [[ !($PROTEIN =~ .*TCR.*) && $PROTEIN =~ .*1bx2.* ]] ; then
cat > stripparm.in << EOF
parm stripped.${PROTEIN}*.parm7
trajin ${PROTEIN}_${Trial}_prod.combine.nc
strip !(@C,CA,N,O) outprefix bb 
trajout bb_${PROTEIN}_${Trial}_1bx2.combine.nc
trajout ${END_DIR}/bb_nontcr_1bx2.nc append
trajout ${END_DIR}/bb_nontcr.nc append
EOF
cpptraj -i stripparm.in > stripparm.log

elif [[ ($PROTEIN =~ .*TCR.*) && $PROTEIN =~ .*5ni9.* ]] ; then
cat > stripparm.in << EOF
parm stripped.${PROTEIN}*.parm7
trajin ${PROTEIN}_${Trial}_prod.combine.nc
strip !(@C,CA,N,O) outprefix bb 
trajout bb_${PROTEIN}_${Trial}_5ni9.combine.nc
trajout ${END_DIR}/backbone_tcr_5ni9.nc append
run
strip !(@C,CA,N,O)|:195-415 outprefix bb 
trajout bbnontcr_${PROTEIN}_${Trial}_5ni9.combine.nc
trajout ${END_DIR}/backbone_tcrbutnontcr_5ni9.nc append
EOF
cpptraj -i stripparm.in > stripparm.log

elif [[ ($PROTEIN =~ .*TCR.*) && $PROTEIN =~ .*1bx2.* ]] ; then
cat > stripparm.in << EOF
parm bb.stripped.${PROTEIN}*.parm7
#trajin ${PROTEIN}_${Trial}_prod.combine.nc
#strip !(@C,CA,N,O) outprefix bb 
#trajout bb_${PROTEIN}_${Trial}_1bx2.combine.nc
trajin bb_${PROTEIN}_${Trial}_1bx2.combine.nc
trajout ${END_DIR}/backbone_tcr_1bx2.nc append
#run
#strip !(@C,CA,N,O)|:195-415 outprefix bb 
#trajout bbnontcr_${PROTEIN}_${Trial}_1bx2.combine.nc
trajin bbnontcr_${PROTEIN}_${Trial}_1bx2.combine.nc
trajout ${END_DIR}/backbone_tcrbutnontcr_1bx2.nc append
run
EOF
cpptraj -i stripparm.in > stripparm.log

else 
# cat > combine_nc.in << EOF
# parm stripped.${PROTEIN}*.parm7
# trajin ${PROTEIN}_prod.0*.nc
# autoimage
# rms ToFirst out ${PROTEIN}_${Trial}_rmsd.dat @N,CA,C,O
# trajout ${PROTEIN}_${Trial}_4x5w.combine.nc
# run
# strip !(@C,CA,N,O) outprefix bb 
# trajout bb_${PROTEIN}_${Trial}_4x5w.combine.nc
# trajout ${END_DIR}/bb_tcr_4x5w.nc append
# run
# EOF
# cpptraj -i combine_nc.in > combine.log
cat > stripparm.in << EOF
parm stripped.${PROTEIN}.parm7
trajin ${PROTEIN}_${Trial}_prod.combine.nc
strip !(@C,CA,N,O) outprefix bb 
trajout bb_${PROTEIN}_${Trial}_4x5w.combine.nc
trajout ${END_DIR}/backbone_tcr_4x5w.nc append
run
strip !(@C,CA,N,O)|:195-415 outprefix bb 
trajout bbnontcr_${PROTEIN}_${Trial}_4x5w.combine.nc
trajout ${END_DIR}/backbone_tcrbutnontcr_4x5w.nc append
run
EOF
cpptraj -i stripparm.in > stripparm.log
fi

DECIDER=`grep "Read 100000 frames and processed 100000 frames." stripparm.log| wc -l `
if ((${DECIDER} >= 1)); then
    echo $PROTEIN $Trial >> ${END_DIR}/succes_concate.txt
else
    echo $dir >> ${END_DIR}/failed_concate.txt
fi
done