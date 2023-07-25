#setenv PATH /dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
#python /dors/meilerlab/home/chey120/mhcii_asyn/scripts/md_analysis/combine_data.py _asyn_rmsf.dat

import sys
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sns
import os
import numpy as np

inp=sys.argv[1]
directory = os.getcwd() 
df_list=[]
for filename in os.listdir(directory):
    # print(filename)
    f = os.path.join(directory, filename)
    if os.path.isfile(f):
        pdb = pd.read_table(filename, sep = "\s+")
        fl = filename.replace(inp, "").split(".00")
        fl1 = fl[0].split("_")
        # print(fl1)
        pdb.insert(1, "Protein", '_'.join(fl1[0:-2]), allow_duplicates = False)
        pdb.insert(1, "Trial", fl1[-1], allow_duplicates = False)
        pdb.insert(1, "MHCII", fl1[-2], allow_duplicates = False)
        df_list.append(pdb)

df = pd.concat(df_list)
df = df.rename(columns={'#Res': 'Res'})
list(df.columns)
out = df[['Res', 'MHCII', 'Trial', 'Protein', 'AtomicFlx']]
out.to_csv('../sum'+inp) 
# print(out.head(100))

# Set working directory
# No equivalent function in Python, as it operates on the file system directly

# Read the data file
# asyn = pd.read_csv("/dors/meilerlab/home/chey120/mhcii_asyn/Data/rmsf/sum_asyn_rmsf.dat")

# Perform data manipulation using pandas
# asyn['Res'] = pd.Series(
#     pd.np.select(
#         [
#             asyn['Res'] > 196,
#             asyn['Res'] <= 196
#         ],
#         [
#             asyn['Res'] - 400,
#             asyn['Res'] - 180
#         ],
#         default=asyn['Res']
#     )
# )
# asyn['Protein'] = pd.Categorical(
#     asyn['Protein'],
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
#     data=asyn, x='Res', y='AtomicFlx', hue='Trial',
#     style='Trial', size='Trial', sizes=[20, 20, 20], legend='full'
# )
# plt.ylim(0, 10)
# plt.xlabel('Res')
# plt.ylabel('AtomicFlx')
# plt.title('RMSF Analysis')
# plt.legend(title='Trial', loc='upper right')
# plt.savefig("AtomicFlx.png")
# plt.show()