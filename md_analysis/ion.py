import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
def ion(file_prefix, column_name, ax):
    dataframes = []
    trials = 5
    for i in range(trials):
        file_name = f"{file_prefix}_trial{i}_{column_name}"
        df = pd.read_csv(file_name, delim_whitespace=True, comment='#', header=None, names=["Distance_(Ang)", "Value"])
        df["Trial"] = i
        dataframes.append(df)

    combined = pd.concat(dataframes)
    grouped_dataset = combined.groupby(['Distance_(Ang)'])['Value']
    average_dataset = grouped_dataset.mean()
    std_dataset = grouped_dataset.std()
    stderr_dataset = std_dataset / np.sqrt(trials)
    x = average_dataset.index
    y = average_dataset.values
    std = stderr_dataset.values

    ax.errorbar(x, y, yerr=std, label=file_prefix, alpha=0.7)
    ax.set_ylim(0, 0.6)
    ax.set_xlabel(column_name.replace(".dat","") + " radial distribution function")
    ax.legend()


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
ion("DNEAY_5ni9", "Na.dat", ax1)
ion("S129_DNEAY_5ni9", "Na.dat", ax1)
ion("DNEAY_3pl6", "Na.dat", ax2)
ion("S129_DNEAY_3pl6", "Na.dat", ax2)
ion("KEGVL_1bx2", "Na.dat", ax3)
ion("Y39_KEGVL_1bx2", "Na.dat", ax3)
ion("KEGVL_1h15", "Na.dat", ax4)
ion("Y39_KEGVL_1h15", "Na.dat", ax4)
plt.savefig(f"ion_stderr.png")
