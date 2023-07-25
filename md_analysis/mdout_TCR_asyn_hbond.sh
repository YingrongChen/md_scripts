#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md/mdout_TCR_process.sh to_process
source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt

# export PATH=/dors/meilerlab/home/chey120/mhcii_asyn/scripts/md:${PATH}
# chmod u+x /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md/mdout_extract.perl
# module load GCC OpenMPI R

#find . -name prod_47381594.log

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
# asyn=":401-415"
# mhcii=":222-400"
# tcr=":1-221"
# cut="5.0"

# # Loop over all proteins in the list
# for dir in `cat ${LIST}`; do
# cd $dir
# Trial=`basename $dir`
# PROTEIN=${dir%*/trial*} #remove /trial0
# PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

# ######
# # extract the data from the output files and plot
# ######
# RESTART=`ls ${PROTEIN}_prod.*.crd | tail -n1 | awk -F. '{print $(NF-1)}'`
# if [ $RESTART == 0100 ]; then

# ######
# # cpptraj: combine NC's to nc Binpos, calculate rmsf of peptide and hydrogen bond
# ######
# cat > combine_nc.in << EOF
# parm stripped.${PROTEIN}*.parm7
# trajin ${PROTEIN}_${Trial}_prod.combine.nc
# autoimage
# rms ToFirst out ${PROTEIN}_${Trial}_rmsd.dat @N,CA,C,O
# hbond donormask ${tcr}@C,CA,N,O acceptormask ${mhcii}@C,CA,N,O avgout ${PROTEIN}_${Trial}_tcrasyn_hbt.dat nointramol
# hbond donormask ${mhcii}@C,CA,N,O acceptormask ${tcr}@C,CA,N,O avgout ${PROTEIN}_${Trial}_tcrasyn_hbt.dat nointramol
# run
# quit
# EOF
# cpptraj -i combine_nc.in &
# fi
# done

for dir in `cat ${LIST}`; do
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/
cp ${PROTEIN}_${Trial}_tcrasyn_hbt.dat /dors/meilerlab/home/chey120/mhcii_asyn/Data/hbond_tcr/asyntcr &
done