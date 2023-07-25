import glob
import pandas as pd
import os

# Get the parent directory name
parent_dir = os.path.basename(os.path.dirname(os.getcwd()))

# Get a list of all files ending with "_score.sc"
file_list = glob.glob("*_score.sc")

# Initialize an empty DataFrame to store the data
data = pd.DataFrame()

# Iterate over each file
for file_name in file_list:
    # Read the content of the current file
    df = pd.read_csv(file_name, delim_whitespace=True, skiprows=1)
    data = pd.concat([data, df], ignore_index=True)

# Sort the DataFrame by the "total_score" column in ascending order
data.sort_values("total_score", ascending=True, inplace=True)

# Write the sorted DataFrame to a new file
output_file = "../../" + parent_dir + ".sc"
data.to_csv(output_file, sep='\t', index=False)
