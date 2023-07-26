import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Define the file pattern
inp = "_secstruct.dat"

# Get a list of all .agr files in the directory
file_list = glob.glob("*"+inp)

# Initialize empty DataFrames
dataframes = []

# Loop through each file
for file_name in file_list:

    df = pd.read_csv(file_name, delim_whitespace=True, skiprows=1, header=None,
                    names=['Residue', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend'])
    fl1 = file_name.replace(inp, "").split("_")
    print(fl1)
    df.insert(1, "Protein", '_'.join(fl1[0:-2]), allow_duplicates = False)
    df.insert(1, "Trial", fl1[-1], allow_duplicates = False)
    df.insert(1, "MHCII", fl1[-2], allow_duplicates = False)
    # Append the DataFrame to the list
    dataframes.append(df)

# Combine the DataFrames into a single DataFrame
combined_df = pd.concat(dataframes, ignore_index=True)

combined_df.to_csv("secstruct.csv", index=False)

residue_67_df = combined_df[combined_df['Residue'] == 67]

# Create subplots
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# Group the DataFrame by MHCII and Protein columns
grouped = residue_67_df.groupby(['MHCII', 'Protein'])

chainA_num = 0
protein_num = 0
# Iterate over the groups and plot stacked bar graph on respective subplot
for (chainA, protein), group in grouped:
    ax = axes[chainA_num, protein_num]  # Calculate subplot index based on MHCII and Protein
    group.drop(columns=['Residue']).set_index('Trial').plot(kind='bar', stacked=True, ax=ax)
    ax.set_title(f'MHCII {chainA} - Protein {protein}')
    ax.legend().set_visible(False)
    # Update subplot index
    protein_num += 1
    if protein_num == 2:
        protein_num = 0
        chainA_num += 1

# Add legend outside the subplots
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right')

plt.tight_layout()
plt.savefig(f"DNEAY_secstruct.png")

# # Group the dataframes by 'InputFile' and 'Label' and sum the 'y' values
# grouped_67 = df_67.groupby(['InputFile', 'Label'])['y'].sum().unstack()
# grouped_65 = df_65.groupby(['InputFile', 'Label'])['y'].sum().unstack()

# # Set the figure size
# fig, (ax_67, ax_65) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# grouped_67.plot(kind='bar', stacked=True, ax=ax_67)
# ax_67.set_xlabel('DRB1*0401 Trajectories')
# ax_67.set_xticklabels(df_67['InputFile'])
# ax_67.legend().remove()
# ax_67.set_ylabel('Secondary Structures')

# grouped_65.plot(kind='bar', stacked=True, ax=ax_65)
# ax_65.set_xlabel('DRB1*0101 Trajectories')
# ax_65.set_xticklabels(df_65['InputFile'])
# ax_65.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
# plt.tight_layout()

# # Save the plot as an image file
# plt.savefig('stacked_bar_charts_combined.png',bbox_inches='tight', dpi=300)
