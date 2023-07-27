import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def rmsf(file_prefix, column_name, ax):
    dataframes = []
    trials = 5
    for i in range(trials):
        file_name = f"{file_prefix}_trial{i}_{column_name}"
        df = pd.read_csv(file_name, delim_whitespace=True, comment='#', header=None, names=["Residue", "Value"])
        df["Trial"] = i
        if "5ni9" in file_prefix:
            df['Value'] = df['Value'].shift(-2)
            df = df.dropna(subset=['Value'])
            df = df.iloc[:-1]
        df = df.iloc[:-1]
        dataframes.append(df)


    combined = pd.concat(dataframes)
    grouped_dataset = combined.groupby(['Residue'])['Value']
    average_dataset = grouped_dataset.mean()
    std_dataset = grouped_dataset.std()
    stderr_dataset = std_dataset / np.sqrt(trials)
    x = average_dataset.index
    y = average_dataset.values
    std = stderr_dataset.values

    ax.errorbar(x, y, yerr=std, label=file_prefix, alpha=0.7)
    ax.set_ylim(0, 6)
    ax.set_xlabel(column_name.replace(".dat",""))
    ax.legend()


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
rmsf("DNEAY_5ni9", "asyn_rmsf.dat", ax1)
rmsf("S129_DNEAY_5ni9", "asyn_rmsf.dat", ax1)
rmsf("DNEAY_3pl6", "asyn_rmsf.dat", ax2)
rmsf("S129_DNEAY_3pl6", "asyn_rmsf.dat", ax2)
rmsf("KEGVL_1bx2", "asyn_rmsf.dat", ax3)
rmsf("Y39_KEGVL_1bx2", "asyn_rmsf.dat", ax3)
rmsf("KEGVL_1h15", "asyn_rmsf.dat", ax4)
rmsf("Y39_KEGVL_1h15", "asyn_rmsf.dat", ax4)
plt.savefig(f"asyn_rmsf_stderr.png")
