# coding=utf-8
#! /usr/bin/env python
import math
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.stats as stats
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

cutoff=3.5
dwell=0
fwd_transitions = []
framesize=0.004 #4fs
sum = 0
input = sys.argv[1]
output = "dwell_"+input
# data = pd.read_csv(input, delim_whitespace=True, on_bad_lines='skip')
# print(data.describe())

with open(input, 'r') as r, open(output, 'a') as w:
    header=r.readline()
    # w.write(input + "\n")
    # w.write("forward:\n")
    for line in r:
        try:
            distance = float(line)
        except ValueError as e:
            continue
        if distance < cutoff:
            dwell = dwell + 1
        if distance >= cutoff and dwell > 0:
            #w.write(str(dwell) + '\n')
            fwd_transitions.append(dwell)
            sum += dwell
            dwell = 0
    # w.write("#frames in forward:" + str(sum) + "\n")

    fwd_transitions = np.array(sorted(fwd_transitions,reverse=True))
    fwd_transition_time = fwd_transitions * framesize
    fwd_transition_prob = []
    for i in range(1, len(fwd_transitions)):
        fwd_transition_prob.append(float(i-1) / float(len(fwd_transitions)))
    fwd_slope, fwd_intercept, fwd_r, fwd_p, fwd_std_err = stats.linregress(fwd_transition_time[1:-1],np.log(fwd_transition_prob[1:]))
    kf_direct = -1.0 * fwd_slope

    # w.write("fwd_slope:" + str(fwd_slope) + "\n")
    # w.write("fwd_intercept:" + str(fwd_intercept) + "\n")
    # w.write("kf_direct:" + str(kf_direct) + "\n")

rev_sum = 0
bwd_transitions = []
with open(input, 'r') as r, open(output, 'a') as w:
    header=r.readline()
    # w.write("reverse:\n")
    for line in r:
        try:
            distance = float(line)
        except ValueError as e:
            continue
        if distance > cutoff:
            dwell = dwell + 1
        if distance <= cutoff and dwell > 0:
            #w.write(str(dwell) + '\n')
            bwd_transitions.append(dwell)
            sum += dwell
            rev_sum += dwell
            dwell = 0
    # w.write("#frames in reverse:" + str(rev_sum) + "\n")
    # w.write("total #frames:" + str(sum) + "\n") 
    # Reorder the transition event time series
    bwd_transitions = np.array(sorted(bwd_transitions,reverse=True))
    bwd_transition_time = bwd_transitions * framesize

    # Compute probability that system escapes energy well at time T or longer
    bwd_transition_prob = []
    for i in range(1, len(bwd_transitions)):
        bwd_transition_prob.append(float(i-1) / float(len(bwd_transitions)))
    bwd_transition_prob[1:]
    bwd_slope, bwd_intercept, bwd_r, bwd_p, bwd_std_err = stats.linregress(bwd_transition_time[1:-1],np.log(bwd_transition_prob[1:]))
    kr_direct = -1.0 * bwd_slope

    # w.write("bwd_slope:" + str(bwd_slope) + "\n")
    # w.write("bwd_intercept:" + str(bwd_intercept) + "\n")
    # w.write("kr_direct:" + str(kr_direct) + "\n")

with open(output, 'a') as w:
    temp = 300 # K ; a standard simulation temperature
    beta = 1.0/(0.001987*temp)
    dG_direct = -1.0 / beta * math.log(kr_direct / kf_direct) #dG = - RTln(keq), ln(n) = 2.303 * log(n), dG = 1.9872036(11)×10−3*temp*2.303*log(n) kcal mol-1
    w.write("k_on - Direct Approach: " + str('{:0.3e}'.format(kf_direct)) + " s-1 \n")
    w.write("k_off - Direct Approach: " + str('{:0.3e}'.format(kr_direct)) + " s-1 \n")
    w.write("dG_rate_constants - Direct Approach: " + str('{:0.3e}'.format(dG_direct)) + " kcal/mol \n")

plt.figure()
plt.subplot(211)
plt.subplots_adjust(hspace = .5, wspace=.3)
plt.plot(fwd_transition_time[1:-1], fwd_transition_prob[1:],marker=".",ms=8,color="blue",linestyle="--",lw=.5)
plt.plot(bwd_transition_time[1:-1], bwd_transition_prob[1:],marker=".",ms=8,color="red",linestyle="--",lw=.5)
plt.xlabel("Time (s)")
plt.ylabel("S(t)")
plt.legend(("Forward Reaction", "Reverse Reaction"),loc="upper right",fontsize=8)
plt.grid(True)
plt.subplot(212)
plt.plot(fwd_transition_time[1:-1], np.log(fwd_transition_prob[1:]),marker=".",ms=8,color="blue",linestyle=":",lw=.5)
plt.plot(bwd_transition_time[1:-1], np.log(bwd_transition_prob[1:]),marker=".",ms=8,color="red",linestyle=":",lw=.5)
plt.plot(fwd_transition_time, fwd_slope*fwd_transition_time + fwd_intercept,color="blue",linestyle="-",lw=1.0)
plt.plot(bwd_transition_time, bwd_slope*bwd_transition_time + bwd_intercept,color="red",linestyle="-",lw=1.0)
plt.xlabel("Time (s)")
plt.ylabel("ln(S(t))")
plt.ylim(-15, 0)
plt.legend(("Forward Reaction", "Reverse Reaction","Forward Reaction Regression","Reverse Reaction Regression"),loc="upper right",fontsize=8)
plt.suptitle(input + " cutoff = " + str(cutoff) + str(" {:0.3e}".format(dG_direct)) + " kcal/mol")
plt.grid(True)
plt.savefig(input[0:-4] + "_cut" + str(cutoff) + "_survival_functions", dpi=300)
plt.close()