for file in `ls *.dat`; do
touch sum_tcr
name=${file%*_trial*}
cat $file | awk -v v=$name '{$(NF+1)=v; print}'>> sum_tcr
done
