import matplotlib.pyplot as plt
import numpy as np
import pyemma
import sys

topfile=('bbtcr.pdb')
trajfile=('bbtcr_small.nc')
feat = pyemma.coordinates.featurizer(topfile)
ser=feat.select("(resSeq 409) and name CA")
tyr=feat.select("(resSeq 405) and name CA")
tyr2=feat.select("(resSeq 413) and name CA")
notser=feat.select("(not resSeq 409) and name CA")
nottyr=feat.select("(not (resSeq 409 or resSeq 405)) and name CA")
nottyr2=feat.select("(not (resSeq 409 or resSeq 405 or resSeq 413)) and name CA")
feat.add_backbone_torsions()
print('number of backbone torsions features ', feat.dimension())
feat.add_distances(notser, indices2=ser)
feat.add_distances(nottyr, indices2=tyr)
feat.add_distances(nottyr2, indices2=tyr2)
print('number of features ', feat.describe())
print('number of features ', feat.dimension())
data = pyemma.coordinates.load(trajfile, features=feat)
tica = pyemma.coordinates.tica(data)
cluster = pyemma.coordinates.cluster_kmeans(tica, k=100, max_iter=50)
its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=3, errors='bayes')
pyemma.plots.plot_implied_timescales(its, outfile="msm_its_ser-notser.png", show_mle=False);
msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=10, dt_traj='10 ps')
pyemma.plots.plot_cktest(msm.cktest(4), units='ps');
plt.savefig("msm_ck_ser-notser.png")

print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))
print(msm.active_set)
print(msm.stationary_distribution)
print('sum of weights = {:f}'.format(msm.pi.sum()))

pyemma.plots.plot_free_energy(*data.T[0:2], legacy=False)
plt.savefig("free_energy_ser-notser.png")
