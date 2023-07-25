import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def graph(input_files, label, data):
    n=0
    for input_file in input_files:
        df = pd.read_csv(input_file, delim_whitespace=True)
        min_values = df.iloc[:, 1:].min(axis=1).values
        df_temp = pd.DataFrame({"Label": label, "Min_Values": min_values, "Trial": n})
        print(n)
        data = pd.concat([data, df_temp], ignore_index=True)
        n=n+1
    return data

input_5ni9 = glob.glob("TCR_S129_DNEAY_5ni9_trial?_188hb397.dat")
input_4x5w = glob.glob("TCR_S129_DNEAY_trial?_409hb203.dat")
input_1bx2 = glob.glob("TCR_S129_DNEAY_1bx2_trial?_187hb396.dat")
labels = ["5ni9", "4x5w", "1bx2"]  # Define labels for each category
data = pd.DataFrame(columns=["Label", "Min_Values", "Trial"])

for inputs, label in zip([input_5ni9, input_4x5w, input_1bx2], labels):
    data = graph(inputs, label, data)

data.to_csv("188hb397.csv", index=False)

input_5ni9 = glob.glob("TCR_S129_DNEAY_5ni9_trial?_159hb258.dat")
input_4x5w = glob.glob("TCR_S129_DNEAY_trial?_380hb64.dat")
input_1bx2 = glob.glob("TCR_S129_DNEAY_1bx2_trial?_158hb257.dat")

labels = ["5ni9", "4x5w", "1bx2"]  # Define labels for each category
data = pd.DataFrame(columns=["Label", "Min_Values", "Trial"])

for inputs, label in zip([input_5ni9, input_4x5w, input_1bx2], labels):
    data = graph(inputs, label, data)

data.to_csv("159hb258.csv", index=False)