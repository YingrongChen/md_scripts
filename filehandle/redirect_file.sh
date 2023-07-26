HIGHRES=/dors/meilerlab/home/chey120/mhcii_asyn/Y39/highres
ANALYSIS=/dors/meilerlab/home/chey120/mhcii_asyn/Y39/analysis
LEN=100 #expected number of lines in file
#rosetta command line
for file in $HIGHRES/*.sc;
do
    s=$(wc -l < $file)
    if (($s > $LEN)); then    
        name=${file##*/} #remove path/
        name=${name%.*} #remove .sc
        cd $ANALYSIS
        if [ ! -d "$ANALYSIS/$name" ]; then
            mkdir $name
            cd $ANALYSIS/$name
            cp $HIGHRES/$name.sc .
            cp $HIGHRES/$name.silent .
            awk 'FNR==2{print $2 "\t" $29 "\t" $30 "\t" $NF}' $name.sc> sort_$name.sc
            sort -n -k2 $name.sc | awk '{print $2 "\t" $29 "\t" $30 "\t" $NF}'| head -n -2 >> sort_$name.sc
        else
            echo $name "exists"
        fi
    else
        echo $file "job unfinished"
    fi
done &