#accre
#bash ../../scripts/md_analysis/secstruct_analysis.sh to_process
module load GCC/6.4.0-2.28 Intel/2017.4.196
module load gnuplot/5.2.2

# Set the desired cbtics values
cbtics_min=0.0
cbtics_step=1.0
cbtics_max=7.0

LIST=`readlink -e $1` # by readlink -e */trial0 > dir_list
for input_file in `cat ${LIST}`; do
        # Get the output file name by replacing the extension of the input file with ".png"
        output_file_name="${input_file%.gnu}.png"

        # Set the output terminal
        output_terminal="pngcairo"

        # Temporary file for storing the modified content
        temp_file=$(mktemp)

        # Add "set term" line at the beginning
        echo "set term $output_terminal" > "$temp_file"

        # Read the input Gnuplot script file, skip "set cbtics" line, and append to temporary file
        grep -v "set cbtics" "$input_file" >> "$temp_file"

        # Append the corrected "set cbtics" line
        echo "set cbtics ${cbtics_min},${cbtics_step},${cbtics_max}" >> "$temp_file"

        # Append the additional lines for output file
        echo "set output \"$output_file_name\"" >> "$temp_file"
        echo "replot" >> "$temp_file"

        # Overwrite the input file with the modified content
        mv "$temp_file" "$input_file"

        echo "Gnuplot script modified successfully. Input file: $input_file"
        echo "Output file: $output_file_name"

        gnuplot $input_file &
done