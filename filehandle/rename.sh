#bash /dors/meilerlab/home/chey120/chainA_chainA/scripts/filehandle/rename.sh

for file in */*
do
# if [[ $file == *_cap_amber* ]]; then
#   mv "$file" "${file/_cap_amber/}"
# fi
# # if [[ $file == TCR_*_trial* ]] && [[ $file != *1bx2* ]]; then
# #   mv "$file" "${file/_trial/_4x5w_trial}"
# # fi
# if [[ $file == *_5ni9_4x5w* ]]; then
#   mv "$file" "${file/_5ni9_4x5w_trial/_5ni9_trial}"
# fi
# if [[ $file == *_peptide_* ]]; then
#   mv "$file" "${file/_peptide_/_chainA_}"
# fi
# if [[ $file == *stripped.* ]]; then
#   mv "$file" "${file/stripped./}"
# fi
if [[ $file == *CR_KEGV* ]]; then
  mv "$file" "${file/TCR_KEGV/TCR_Y39_KEGV}"
fi
done