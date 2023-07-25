import math
import sys
import scipy.stats as stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

cutoff = 4
framesize = 0.004
input_files = sys.argv[1:]  # List of input files

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))
plt.subplots_adjust(hspace=0.5)

output_file = "dwell_time_dG.csv"
results_df = pd.DataFrame(columns=["Input File", "k_on", "k_off", "dG_rate_constants"])

for input_file in input_files:
    df = pd.read_csv(input_file, delim_whitespace=True)
    min_values = df.iloc[:, 1:].min(axis=1).values
    dwell = 0
    fwd_transitions = []
    sum = 0    
    for distance in min_values:
        if math.isnan(distance):
            continue
        if distance < cutoff:
            dwell += 1
        if distance >= cutoff and dwell > 0:
            fwd_transitions.append(dwell)
            sum += dwell
            dwell = 0
    fwd_transitions = np.array(sorted(fwd_transitions, reverse=True))
    fwd_transition_time = fwd_transitions * framesize
    fwd_transition_prob = []
    for i in range(1, len(fwd_transitions)):
        fwd_transition_prob.append(float(i-1) / float(len(fwd_transitions)))
    fwd_slope, fwd_intercept, fwd_r, fwd_p, fwd_std_err = stats.linregress(fwd_transition_time[1:-1],np.log(fwd_transition_prob[1:]))
    kf_direct = -1.0 * fwd_slope
    
    rev_sum = 0
    bwd_transitions = []
    for distance in min_values:
        if math.isnan(distance):
            continue
        if distance > cutoff:
            dwell += 1
        if distance <= cutoff and dwell > 0:
            bwd_transitions.append(dwell)
            sum += dwell
            rev_sum += dwell
            dwell = 0
    
    bwd_transitions = np.array(sorted(bwd_transitions,reverse=True))
    bwd_transition_time = bwd_transitions * framesize

    bwd_transition_prob = []
    for i in range(1, len(bwd_transitions)):
        bwd_transition_prob.append(float(i-1) / float(len(bwd_transitions)))
    bwd_transition_prob[1:]
    bwd_slope, bwd_intercept, bwd_r, bwd_p, bwd_std_err = stats.linregress(bwd_transition_time[1:-1],np.log(bwd_transition_prob[1:]))
    kr_direct = -1.0 * bwd_slope
    
    ax1.plot(bwd_transition_time[1:-1], bwd_transition_prob[1:],marker=".",ms=1,color="red",alpha = 0.8,linestyle="--",lw=.5)
    
    ax2.plot(fwd_transition_time[1:-1], np.log(fwd_transition_prob[1:]),marker=".",ms=1,color="blue",alpha = 0.8, linestyle=":",lw=.5)
    ax2.plot(bwd_transition_time[1:-1], np.log(bwd_transition_prob[1:]),marker=".",ms=1,color="red",alpha = 0.8,linestyle=":",lw=.5)
    ax2.plot(fwd_transition_time, fwd_slope*fwd_transition_time + fwd_intercept,color="blue",alpha = 0.8,linestyle="-",lw=1.0,label=str(i))
    ax2.plot(bwd_transition_time, bwd_slope*bwd_transition_time + bwd_intercept,color="red",alpha = 0.8,linestyle="-",lw=1.0,label=str(i))
    
    temp = 300 # K ; a standard simulation temperature
    beta = 1.0/(0.001987*temp)
    dG_direct = -1.0 / beta * math.log(kr_direct / kf_direct) #dG = - RTln(keq), ln(n) = 2.303 * log(n), dG = 1.9872036(11)×10−3*temp*2.303*log(n) kcal mol-1
    result = {
        "Input File": input_file,
        "k_on": kf_direct,
        "k_off": kr_direct,
        "dG_rate_constants": dG_direct
    }

    # Append the result to the DataFrame
    results_df = results_df.append(result, ignore_index=True)

ax1.set_ylabel("S(t)")
ax1.legend(("Forward Reaction", "Reverse Reaction"),loc="upper right",fontsize=8)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("ln(S(t))")
ax2.set_ylim(-15, 0)
ax2.set_xlim(0,30)
ax2.legend(("Forward Reaction", "Reverse Reaction","Forward Reaction Regression","Reverse Reaction Regression"),loc="upper right",fontsize=8)
output_filename = os.path.splitext(input_files[0])[0] + ".png"
plt.savefig(output_filename, dpi=300)
plt.close()

file_exists = os.path.isfile(output_file)
with open(output_file, 'a') as f:
    results_df.to_csv(f, header=not file_exists, index=False)