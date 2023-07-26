import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

inp="_backbond_hbt.dat"
directory = os.getcwd() 
df_list=[]
file_names = glob.glob("*" + inp)
for filename in file_names:
    pdb = pd.read_table(filename, sep = "\s+")
    pdb['Frac'] = pd.to_numeric(pdb['Frac'], errors='coerce')
    pdb = pdb[pdb['Frac'] > 0.01]
    pdb = pdb.rename(columns={"#Acceptor": "Acceptor"})
    pdb['acceptorid'] = pd.to_numeric(pdb['Acceptor'].str.split('_').str[1].str.split('@').str[0], errors='coerce')
    pdb['donorid'] = pd.to_numeric(pdb['Donor'].str.split('_').str[1].str.split('@').str[0], errors='coerce')
    pdb['peptideid'] = np.maximum(pdb['acceptorid'], pdb['donorid'])
    pdb['Hbond'] = np.where(pdb['acceptorid'] > pdb['donorid'],
                            pdb['Acceptor'] + "-" + pdb['DonorH'],
                            pdb['DonorH'] + "-" + pdb['Acceptor'])
    fl1 = filename.replace(inp, "").split("_")
    # print(fl1)
    pdb.insert(1, "Protein", '_'.join(fl1[0:-2]), allow_duplicates = False)
    pdb.insert(1, "Trial", fl1[-1], allow_duplicates = False)
    pdb.insert(1, "MHCII", fl1[-2], allow_duplicates = False)
    df_list.append(pdb)

df = pd.concat(df_list)
list(df.columns)
df.to_csv('../sum'+inp, index=False) 
# print(out.head(100))

# Set working directory
# No equivalent function in Python, as it operates on the file system directly

# Read the data file
# chainA = pd.read_csv("/dors/meilerlab/home/chey120/chainA_chainA/Data/rmsf/sum_chainA_rmsf.dat")

# Perform data manipulation using pandas
# chainA['Res'] = pd.Series(
#     pd.np.select(
#         [
#             chainA['Res'] > 196,
#             chainA['Res'] <= 196
#         ],
#         [
#             chainA['Res'] - 400,
#             chainA['Res'] - 180
#         ],
#         default=chainA['Res']
#     )
# )
# chainA['Protein'] = pd.Categorical(
#     chainA['Protein'],
#     categories=[
#         "DNEAY", "Y125_DNEAY", "S129_DNEAY", "TCR_DNEAY", "TCR_Y125_DNEAY",
#         "TCR_S129_DNEAY", "KEGVL", "Y39_KEGVL", "S42_KEGVL", "TCR_KEGVL",
#         "TCR_Y39_KEGVL", "TCR_S42_KEGVL", "AEAAG", "KEGVV"
#     ]
# )

# # Generate the plot using seaborn
# plt.figure(figsize=(10, 5))
# sns.set(style="ticks")
# sns.set_palette("husl")
# sns.scatterplot(
#     data=chainA, x='Res', y='AtomicFlx', hue='Trial',
#     style='Trial', size='Trial', sizes=[20, 20, 20], legend='full'
# )
# plt.ylim(0, 10)
# plt.xlabel('Res')
# plt.ylabel('AtomicFlx')
# plt.title('RMSF Analysis')
# plt.legend(title='Trial', loc='upper right')
# plt.savefig("AtomicFlx.png")
# plt.show()

#setenv PATH /dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
#python /combine_data.py _chainA_rmsf.dat