# bash /dors/meilerlab/home/chey120/mhcii_asyn/scripts/filehandle/dat_combine.sh ../nativecontacts_tcr_asyn.dat
# for file in `ls *.dat`; do
# TRIAL=`basename $file _*_*.dat`
# PROTEIN=${file%*_trial*} #remove /trial0
# cat $file | while read line; do echo -e ${line}' \t '$PROTEIN' \t '$TRIAL >> $filename; done
# # sed -e "s/$/\t$PROTEIN/\t$TRIAL" $file >> $filename
# done &
#sed 's/Y125_/Y125/g;s/Y39_/Y39/g;s/S129_/S129/g;s/S42_/S42/g;s/TCR_/TCR/g' $file -i

for file in `ls */*.pdb`; do
echo $file
sed 's/Y125_/Y125/g' $file -i
sed 's/Y39_/Y39/g' $file -i
sed 's/S129_/S129/g' $file -i
sed 's/S42_/S42/g' $file -i
sed 's/TCR_/TCR/g' $file -i
done &