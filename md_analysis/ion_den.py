import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Get the current working directory
current_directory = os.getcwd()

# Read the file containing directories
with open("to_process", "r") as file:
    directories = file.read().splitlines()

# Initialize an empty DataFrame for the big dataset
big_dataset = pd.DataFrame()

# Initialize an empty list for storing unique filenames
unique_filenames = []

# Loop through each directory
for directory in directories:
    try:
        # Change into the directory
        os.chdir(directory)
        
        # Check if the file exists
        if not os.path.exists("200mM-Na.dat"):
            print(f"File 200mM-Na.dat not found in {directory}. Skipping...")
            os.chdir("..")
            continue
        
        # Extract the second-to-last directory name
        filename = os.path.basename(os.path.dirname(directory))
        
        # Combine filenames ending with "_4x5w" and "_5ni9" into a single group
        if filename.endswith("_4x5w") or filename.endswith("_5ni9"):
            filename = filename[:-5]  # Remove the "_4x5w" or "_5ni9" suffix
        
        # Read the file using pandas
        df = pd.read_csv("200mM-Na.dat", delim_whitespace=True, comment='#', header=None, names=["Distance_(Ang)", "Value"])
        
        # Add the filename as a new column
        df['filename'] = filename
        
        # Append the current data to the big dataset
        big_dataset = big_dataset.append(df, ignore_index=True)
        
        # Append the filename to the list of unique filenames
        if filename not in unique_filenames:
            unique_filenames.append(filename)
    
    except FileNotFoundError:
        print(f"File not found error occurred in {directory}. Skipping...")
        os.chdir("..")
        continue

# Plot the graph for each unique filename
grouped_dataset = big_dataset.groupby(['filename', 'Distance_(Ang)'])['Value']
average_dataset = grouped_dataset.mean()
std_dataset = grouped_dataset.std()

# Plot the graph for each filename
# Plot the graph for each filename
for filename in unique_filenames:
    subset = average_dataset.loc[filename]
    x = subset.index.get_level_values('Distance_(Ang)')
    y = subset.values
    std = std_dataset.loc[filename].values
    plt.errorbar(x, y, yerr=std, label=filename, alpha=0.7)

plt.xlabel("Distance_(Ang)")
plt.ylabel("Value")
plt.title("Na+ density around serine or phospho-S at :189")

# Set the y-axis limits
plt.ylim(0, 0.7)
plt.xlim(0, 10)
# Display the legend
plt.legend()

# Save the figure to the original directory
os.chdir(current_directory)
save_directory = os.path.basename(os.path.dirname(os.path.dirname(directories[-2])))
print(save_directory)
print(current_directory)

plt.savefig(f"ionden_std.png")
plt.clf()

n_samples = grouped_dataset.size().groupby('filename').mean()  # Average number of samples per group
stderr_dataset = std_dataset / np.sqrt(n_samples)

# Plot the graph for each filename
for filename in unique_filenames:
    subset = average_dataset.loc[filename]
    x = subset.index.get_level_values('Distance_(Ang)')
    y = subset.values
    stderr = stderr_dataset.loc[filename].values
    plt.errorbar(x, y, yerr=stderr, label=filename, alpha=0.7)
plt.xlabel("Distance_(Ang)")
plt.ylabel("Value")
plt.title("Na+ density around serine or phospho-S at :189")

# Set the y-axis limits
plt.ylim(0, 0.7)
plt.xlim(0, 10)
# Display the legend
plt.legend()
plt.savefig(f"ionden_stderr.png")


