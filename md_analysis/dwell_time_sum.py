#python md/dwell_time_sum.py *_hbond_188.dat DNEAY_4x5w_sum_hbond_188_sum.dat
import math
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt

cutoff=3.5000
framesize=0.004 #4fs
inputs = sys.argv[1:-1]
output = sys.argv[-1]
fwd_transition_times = []
fwd_transition_probs = []
bwd_transition_times = []
bwd_transition_probs = []
fwd_slopes = []
bwd_slopes = []
fwd_intercepts = []
bwd_intercepts = []

for input in inputs:
    dwell=0
    fwd_transitions = []
    sum = 0
    with open(input, 'r') as r, open(output, 'a') as w:
        header=r.readline()
        for line in r:

            distance = float(line)
            if distance < cutoff:
                dwell = dwell + 1
            if distance >= cutoff and dwell > 0:
                fwd_transitions.append(dwell)
                sum += dwell
                dwell = 0
        fwd_transitions = np.array(sorted(fwd_transitions,reverse=True))
        fwd_transition_time = fwd_transitions * framesize
        fwd_transition_prob = []
        for i in range(1, len(fwd_transitions)):
            fwd_transition_prob.append(float(i-1) / float(len(fwd_transitions)))
        fwd_slope, fwd_intercept, fwd_r, fwd_p, fwd_std_err = stats.linregress(fwd_transition_time[1:-1],np.log(fwd_transition_prob[1:]))
        kf_direct = -1.0 * fwd_slope

    rev_sum = 0
    bwd_transitions = []
    with open(input, 'r') as r, open(output, 'a') as w:
        header=r.readline()
        # w.write("reverse:\n")
        for line in r:
            distance = float(line)
            if distance > cutoff:
                dwell = dwell + 1
            if distance <= cutoff and dwell > 0:
                bwd_transitions.append(dwell)
                sum += dwell
                rev_sum += dwell
                dwell = 0

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

        temp = 300 # K ; a standard simulation temperature
        beta = 1.0/(0.001987*temp)
        dG_direct = -1.0 / beta * math.log(kr_direct / kf_direct) #dG = - RTln(keq), ln(n) = 2.303 * log(n), dG = 1.9872036(11)×10−3*temp*2.303*log(n) kcal mol-1
        w.write(header + "  " + input + "\n")
        w.write("k_on: " + str('{:0.3}'.format(kf_direct)) + " s-1 \n")
        w.write("k_off: " + str('{:0.3}'.format(kr_direct)) + " s-1 \n")
        w.write("dG_rate_constants: " + str('{:0.3}'.format(dG_direct)) + " kcal/mol \n \n")

        fwd_transition_times.append(fwd_transition_time)
        fwd_transition_probs.append(fwd_transition_prob)
        bwd_transition_times.append(bwd_transition_time)
        bwd_transition_probs.append(bwd_transition_prob)
        fwd_slopes.append(fwd_slope)
        bwd_slopes.append(bwd_slope)
        fwd_intercepts.append(fwd_intercept)
        bwd_intercepts.append(bwd_intercept)

plt.figure()
for i in range(0, len(fwd_transition_times)):
    plt.subplot(211)
    plt.subplots_adjust(hspace = .5, wspace=.3)
    plt.plot(fwd_transition_times[i][1:-1], fwd_transition_probs[i][1:],marker=".",ms=1,color="blue",alpha = 0.8,linestyle="--",lw=.5)
    plt.plot(bwd_transition_times[i][1:-1], bwd_transition_probs[i][1:],marker=".",ms=1,color="red",alpha = 0.8,linestyle="--",lw=.5)
    plt.xlabel("Time (s)")
    plt.ylabel("S(t)")
    plt.legend(("Forward Reaction", "Reverse Reaction"),loc="upper right",fontsize=8)
    plt.grid(True)
    plt.subplot(212)
    plt.plot(fwd_transition_times[i][1:-1], np.log(fwd_transition_probs[i][1:]),marker=".",ms=1,color="blue",alpha = 0.8, linestyle=":",lw=.5)
    plt.plot(bwd_transition_times[i][1:-1], np.log(bwd_transition_probs[i][1:]),marker=".",ms=1,color="red",alpha = 0.8,linestyle=":",lw=.5)
    plt.plot(fwd_transition_times[i], fwd_slopes[i]*fwd_transition_times[i] + fwd_intercepts[i],color="blue",alpha = 0.8,linestyle="-",lw=1.0,label=str(i))
    plt.plot(bwd_transition_times[i], bwd_slopes[i]*bwd_transition_times[i] + bwd_intercepts[i],color="red",alpha = 0.8,linestyle="-",lw=1.0,label=str(i))
    plt.xlabel("Time (s)")
    plt.ylabel("ln(S(t))")
    plt.legend(("Forward Reaction", "Reverse Reaction","Forward Reaction Regression","Reverse Reaction Regression"),loc="upper right",fontsize=8)
    plt.grid(True)
plt.title(output[0:-4])
plt.savefig(output[0:-4] + "_survival_functions", dpi=300)
plt.close()