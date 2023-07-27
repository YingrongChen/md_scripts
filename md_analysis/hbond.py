import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob

inp="backbond_hbt"
def hbond(file_prefix, column_name, ax):
    dataframes = []
    trials = 5
    for i in range(trials):
        file_name = "{}_trial{}_{}.dat".format(file_prefix, i, column_name)
        pdb = pd.read_table(file_name, sep = "\s+")
        pdb['Frac'] = pd.to_numeric(pdb['Frac'], errors='coerce')
        pdb = pdb[pdb['Frac'] > 0.01]
        pdb = pdb.rename(columns={"#Acceptor": "Acceptor"})
        pdb['acceptorid'] = pd.to_numeric(pdb['Acceptor'].str.split('_').str[1].str.split('@').str[0], errors='coerce')
        pdb['donorid'] = pd.to_numeric(pdb['Donor'].str.split('_').str[1].str.split('@').str[0], errors='coerce')
        pdb['peptideid'] = np.maximum(pdb['acceptorid'], pdb['donorid'])
        if "5ni9" in file_prefix:
            pdb['peptideid'] = pdb['peptideid'] - 2
        pdb = pdb[pdb['peptideid'] >= 180]
        pdb['Hbond'] = np.where(pdb['acceptorid'] > pdb['donorid'],
                                f"{pdb['acceptorid']}_{pdb['donorid']}",
                                f"{pdb['donorid']}_{pdb['acceptorid']}")
        dataframes.append(pdb)

    combined = pd.concat(dataframes)
    hbond_counts = combined['Hbond'].value_counts()
    hbond_frequent = hbond_counts[hbond_counts > 3].index
    combined = combined[combined['Hbond'].isin(hbond_frequent)]
    hbond_stats = combined.groupby(['Hbond', 'peptideid'])['Frac'].agg(['mean', 'sem']).reset_index()
    legend_label = file_prefix.replace("_", " ")
    ax.errorbar(hbond_stats['peptideid'], hbond_stats['mean'], yerr=hbond_stats['sem'], fmt='o', label=legend_label)
    ax.legend(loc='upper left')
    ax.set_ylim(0,1)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
hbond("DNEAY_5ni9", inp, ax1)
hbond("S129_DNEAY_5ni9", inp, ax1)
hbond("DNEAY_3pl6", inp, ax2)
hbond("S129_DNEAY_3pl6", inp, ax2)
hbond("KEGVL_1bx2", inp, ax3)
hbond("Y39_KEGVL_1bx2", inp, ax3)
hbond("KEGVL_1h15", inp, ax4)
hbond("Y39_KEGVL_1h15", inp, ax4)
plt.savefig(f"{inp}.png")  

def AvgDist(file_prefix, column_name, ax):
    dataframes = []
    trials = 5
    for i in range(trials):
        file_name = "{}_trial{}_{}.dat".format(file_prefix, i, column_name)
        pdb = pd.read_table(file_name, sep = "\s+")
        pdb['AvgDist'] = pd.to_numeric(pdb['AvgDist'], errors='coerce')
        pdb = pdb[pdb['AvgDist'] > 0.01]
        pdb = pdb.rename(columns={"#Acceptor": "Acceptor"})
        pdb['acceptorid'] = pd.to_numeric(pdb['Acceptor'].str.split('_').str[1].str.split('@').str[0], errors='coerce')
        pdb['donorid'] = pd.to_numeric(pdb['Donor'].str.split('_').str[1].str.split('@').str[0], errors='coerce')
        pdb['peptideid'] = np.maximum(pdb['acceptorid'], pdb['donorid'])
        if "5ni9" in file_prefix:
            pdb['peptideid'] = pdb['peptideid'] - 2
        pdb = pdb[pdb['peptideid'] >= 180]
        pdb['Hbond'] = np.where(pdb['acceptorid'] > pdb['donorid'],
                                f"{pdb['acceptorid']}_{pdb['donorid']}",
                                f"{pdb['donorid']}_{pdb['acceptorid']}")
        dataframes.append(pdb)

    combined = pd.concat(dataframes)
    hbond_counts = combined['Hbond'].value_counts()
    hbond_frequent = hbond_counts[hbond_counts > 3].index
    combined = combined[combined['Hbond'].isin(hbond_frequent)]
    hbond_stats = combined.groupby(['Hbond', 'peptideid'])['AvgDist'].agg(['mean', 'sem']).reset_index()
    legend_label = file_prefix.replace("_", " ")
    ax.errorbar(hbond_stats['peptideid'], hbond_stats['mean'], yerr=hbond_stats['sem'], fmt='o', label=legend_label)
    ax.legend(loc='upper left')
    ax.set_ylim(2.5,3)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
AvgDist("DNEAY_5ni9", inp, ax1)
AvgDist("S129_DNEAY_5ni9", inp, ax1)
AvgDist("DNEAY_3pl6", inp, ax2)
AvgDist("S129_DNEAY_3pl6", inp, ax2)
AvgDist("KEGVL_1bx2", inp, ax3)
AvgDist("Y39_KEGVL_1bx2", inp, ax3)
AvgDist("KEGVL_1h15", inp, ax4)
AvgDist("Y39_KEGVL_1h15", inp, ax4)
plt.savefig(f"{inp}_distance.png")  