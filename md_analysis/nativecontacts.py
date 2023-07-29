import sys
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sns
import os
import numpy as np
import glob


inp="_nativecontacts.dat"
directory = os.getcwd() 
df_list=[]
file_names = glob.glob("*" + inp)
for filename in file_names:
    pdb = pd.read_table(filename, sep = "\s+", comment="#", names=['Num', 'Contact', 'Nframes', "Frac.","Avg","Stdev"])
    pdb[['res1', 'res2']] = pdb['Contact'].str.split('_', expand=True)
    pdb['res1'] = pdb['res1'].str.replace(':', '').str.split('@').str[0]
    pdb['res2'] = pdb['res2'].str.replace(':', '').str.split('@').str[0]

    fl1 = filename.replace(inp, "").split("_")
    pdb.insert(1, "Protein", '_'.join(fl1[0:-2]), allow_duplicates = False)
    pdb.insert(1, "Trial", fl1[-1], allow_duplicates = False)
    pdb.insert(1, "MHCII", fl1[-2], allow_duplicates = False)
    df_list.append(pdb)

df = pd.concat(df_list)
grouped_df = df.groupby(['Protein', 'Trial', 'MHCII', 'res1'])['Frac.'].sum().reset_index()
list(df.columns)
df.to_csv('../sum'+inp, index=False) 
grouped_df.to_csv('../group'+inp, index=False) 