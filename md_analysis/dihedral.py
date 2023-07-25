import pandas as pd
import matplotlib.pyplot as plt

def plot_density(file_prefix, column_name, ax):
    dataframes = []
    for i in range(5):
        file_name = f"{file_prefix}_trial{i}_{column_name}"
        file_path = f"dihedral/{file_name}.dat"
        df = pd.read_table(file_path, sep="\s+", comment="#", names=["Frame", "chi"])
        df["Trial"] = i
        dataframes.append(df)

    combined = pd.concat(dataframes)
    print(combined)

    for trial in combined["Trial"].unique():
        df = combined[combined["Trial"] == trial]
        pd.to_numeric(df['chi']).plot.kde(ax=ax, label=f"Trial {trial}", xlim=(-180, 180), ylim=(0.0, 0.08))
        plt.margins(y=0.2)
    ax.set_xlabel(file_prefix + " " + column_name)
    ax.legend()

# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
# plot_density("KEGVL_1bx2", "tyrchi1", ax1)
# plot_density("Y39_KEGVL_1bx2", "tyrchi1", ax2)
# plot_density("KEGVL_1bx2", "tyrchi2", ax3)
# plot_density("Y39_KEGVL_1bx2", "tyrchi2", ax4)
# plt.tight_layout()
# plt.savefig(f"KEGVL_1bx2_dihedral.png")

# plt.clf()

# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))
# plot_density("KEGVL_1h15", "tyrchi1", ax1)
# plot_density("Y39_KEGVL_1h15", "tyrchi1", ax2)
# plot_density("KEGVL_1h15", "tyrchi2", ax3)
# plot_density("Y39_KEGVL_1h15", "tyrchi2", ax4)
# plt.tight_layout()
# plt.savefig(f"KEGVL_1h15_dihedral.png")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
plot_density("DNEAY_5ni9", "tyrchi1", ax1)
plot_density("S129_DNEAY_5ni9", "tyrchi1", ax2)
plt.tight_layout()
plt.savefig(f"DNEAY_5ni9_dihedral.png")

# # Read and combine KEGVL_1bx2 data frames
# kegvl_dataframes = []
# for i in range(5):
#     file_name = f"KEGVL_1bx2_trial{i}_tyrchi1"
#     file_path = f"dihedral/{file_name}.dat"
#     df = pd.read_table(file_path, sep = "\s+", comment="#", names=["Frame", "chi1"])
#     df["Trial"] = i
#     kegvl_dataframes.append(df)

# kegvl_combined = pd.concat(kegvl_dataframes)
# print(kegvl_combined)
# # Read and combine Y39_KEGVL_1bx2 data frames
# y39_dataframes = []
# for i in range(5):
#     file_name = f"Y39_KEGVL_1bx2_trial{i}_tyrchi1"
#     file_path = f"dihedral/{file_name}.dat"
#     df = pd.read_table(file_path, sep = "\s+", comment="#", names=["Frame", "chi1"])
#     df["Trial"] = i
#     y39_dataframes.append(df)

# y39_combined = pd.concat(y39_dataframes)
# print(y39_combined)

# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))
# for trial in kegvl_combined["Trial"].unique():
#     df = kegvl_combined[kegvl_combined["Trial"] == trial]
#     # plt.hist(df["chi1"], label=f"Trial {trial}", bins=180)
#     pd.to_numeric(df['chi1']).plot.kde(ax=ax1, label=f"Trial {trial}", xlim=(-180, 180))

# for trial in y39_combined["Trial"].unique():
#     df = y39_combined[y39_combined["Trial"] == trial]
#     # plt.hist(df["chi1"], label=f"Trial {trial}", bins=180)
#     pd.to_numeric(df['chi1']).plot.kde(ax=ax2, label=f"Trial {trial}", xlim=(-180, 180))

# plt.savefig("Y39_1bx2_dihedral.png")
