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