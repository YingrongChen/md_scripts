#bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md/hbonddistance_DNEAY.sh dir_list 
#source /dors/meilerlab/home/brownbp1/PowerCoding_10-24-2020/MD_Simulations/environment_control/amber18_environment.bash.txt
# module load GCC OpenMPI R
forw=/dors/meilerlab/home/chey120/mhcii_asyn/Data/DNEAY_forw.txt
rev=/dors/meilerlab/home/chey120/mhcii_asyn/Data/DNEAY_rev.txt
LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list

# Loop over all proteins in the list
for dir in `cat ${LIST}`; do
cd $dir
Trial=`basename $dir`
PROTEIN=${dir%*/trial*} #remove /trial0
PROTEIN=${PROTEIN#*prod/*} #remove things before prod/

# cat > distance_hbond.in << EOF
# parm stripped.${PROTEIN}*.parm7 
# trajin ${PROTEIN}_prod.combine.binpos
# EOF

# if [[ $PROTEIN =~ "rev" ]]; then
# while read -r line; do
# hbond=`echo ${line}| tr -d '[:space:]' | tr -d ':'`
# cat >> distance_hbond.in << EOF
# distance ${line} out distance_hbond.dat
# EOF
# done < $rev

# else
# while read -r line; do
# hbond=`echo ${line}| tr -d '[:space:]'`
# cat >> distance_hbond.in << EOF
# distance ${line} out distance_hbond.dat
# EOF
# done < $forw
# fi
# cpptraj -i distance_hbond.in

#post-process
for dat in hb_*.dat; do
hbond=${dat%*.dat}
hbond=${hbond#hb_*} 
head=`head -n 1 distance_hbond.dat | sed "s/#//g"`
if [[ $PROTEIN =~ "rev" ]]; then
for i in `seq 1 1 14`; do
hbond=`head -n $i $rev | tail -1 | tr -d '[:space:]' | tr -d ':'`
i=`echo $i | awk '{printf( "%05d", $0)}'`
head="${head/Dis_${i}/"$hbond"}"
done

else
for i in `seq 1 1 18`; do
hbond=`head -n $i $forw | tail -1 | tr -d '[:space:]' | tr -d ':'`
i=`echo $i | awk '{printf( "%05d", $0)}'`
head="${head/Dis_${i}/"$hbond"}"
done
fi
done
sed "1s/.*/$head/" distance_hbond.dat > temp_distance_hbond.dat
# Rscript /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md/hbonddistance_plot.R temp_distance_hbond.dat $PROTEIN $Trial
# if [[ $Trial =~ "trial0" ]]; then
# cat temp_distance_hbond.dat > ../../${PROTEIN}_hbond.dat
# else
# tail -n+2 distance_hbond.dat >> ../../${PROTEIN}_hbond.dat
# fi
done &
