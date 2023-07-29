import pandas as pd
import matplotlib.pyplot as plt
def width(file_prefix, column_name, ax):
    dataframes = []
    for i in range(5):
        file_name = f"{file_prefix}_trial{i}_{column_name}.dat"
        df = pd.read_table(file_name, sep="\s+", comment="#", names=["Frame", "Dis"])
        df["Trial"] = i
        dataframes.append(df)

    combined = pd.concat(dataframes)
    print(combined)

    for trial in combined["Trial"].unique():
        df = combined[combined["Trial"] == trial]
        pd.to_numeric(df['Dis']).plot.kde(ax=ax, label=f"Trial {trial}")
        plt.margins(y=0.2)
    ax.set_xlabel(file_prefix + " " + column_name)
    ax.legend()
    ax.set_xlim(15,30)
    ax.set_ylim(0,1)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
width("KEGVL_1bx2", "mhciiwidth", ax1)
width("Y39_KEGVL_1bx2", "mhciiwidth", ax2)
width("KEGVL_1h15", "mhciiwidth", ax3)
width("Y39_KEGVL_1h15", "mhciiwidth", ax4)
plt.tight_layout()
plt.savefig(f"KEGVL_width.png")

# plt.clf()

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
width("DNEAY_5ni9", "mhciiwidth", ax1)
width("S129_DNEAY_5ni9", "mhciiwidth", ax2)
width("DNEAY_3pl6", "mhciiwidth", ax3)
width("S129_DNEAY_3pl6", "mhciiwidth", ax4)
plt.tight_layout()
plt.savefig(f"DNEAY_width.png")