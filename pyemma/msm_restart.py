
#export PATH=/dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
matplotlib.rcParams.update({'font.size': 12})
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import warnings


topfile=sys.argv[1]
trajfile=sys.argv[2]
nstates = 4

def plot_sampled_function(xall, yall, zall, ax=None, nbins=100, nlevels=20, cmap=plt.cm.bwr, cbar=True, cbar_label=None):
    # histogram data
    xmin = np.min(xall)
    xmax = np.max(xall)
    dx = (xmax - xmin) / float(nbins)
    ymin = np.min(yall)
    ymax = np.max(yall)
    dy = (ymax - ymin) / float(nbins)
    # bin data
    eps = x
    xbins = np.linspace(xmin - 0.5*dx, xmax + 0.5*dx, num=nbins)
    ybins = np.linspace(ymin - 0.5*dy, ymax + 0.5*dy, num=nbins)
    xI = np.digitize(xall, xbins)
    yI = np.digitize(yall, ybins)
    # result
    z = np.zeros((nbins, nbins))
    N = np.zeros((nbins, nbins))
    # average over bins
    for t in range(len(xall)):
        z[xI[t], yI[t]] += zall[t]
        N[xI[t], yI[t]] += 1.0

    with warnings.catch_warnings() as cm:
        warnings.simplefilter('ignore')
        z /= N
    # do a contour plot
    extent = [xmin, xmax, ymin, ymax]
    if ax is None:
        ax = plt.gca()
    ax.contourf(z.T, 100, extent=extent, cmap=cmap)
    if cbar:
        cbar = plt.colorbar()
        if cbar_label is not None:
            cbar.ax.set_ylabel(cbar_label)

    return ax

def plot_sampled_density(xall, yall, zall, ax=None, nbins=100, cmap=plt.cm.Blues, cbar=True, cbar_label=None):
    return plot_sampled_function(xall, yall, zall, ax=ax, nbins=nbins, cmap=cmap, cbar=cbar, cbar_label=cbar_label)

feat = pyemma.coordinates.featurizer(topfile)
ser=feat.select("(resSeq 188) and name CA")
tyr=feat.select("(resSeq 184) and name CA")
tyr2=feat.select("(resSeq 192) and name CA")
# nottail="(resSeq 4 to 81 and 88 to 175)"
notser=feat.select("(not resSeq 188) and name CA")
nottyr=feat.select("(not (resSeq 188 or resSeq 184)) and name CA")
nottyr2=feat.select("(not (resSeq 188 or resSeq 184 or resSeq 192)) and name CA")

ahelix=feat.select("(resSeq 51 to 77) and name CA")
bhelix=feat.select("(resSeq 139 to 172) and name CA")
peptide=feat.select("resSeq 181 to 194 and name CA")

n = len(ahelix)
m = len(peptide)
matrix = np.zeros((3*n + 6*m,2), dtype=int)
for j in range(0,3):
    for i in range(0,n):
        matrix[i+j*n,:] = [bhelix[i+j], ahelix[n-1-i]]
for j in range(0,3):
    for i in range(0,m):
        matrix[i+3*n+j*m,:] = [peptide[i], ahelix[n-1-i-j*2]]
for j in range(0,3):
    for i in range(0,m):
        matrix[i+3*n+(3+j)*m,:] = [peptide[i], bhelix[n-1-i-j*2]]
feat.add_backbone_torsions()
feat.add_distances(notser, indices2=ser)
feat.add_distances(nottyr, indices2=tyr)
feat.add_distances(nottyr2, indices2=tyr2)
feat.add_distances(matrix)
print('number of features ', feat.dimension())

data = pyemma.coordinates.source(trajfile, features=feat)
tica = pyemma.load('tica.pyemma', model_name='tica')
Y = tica.get_output() # get tica coordinates
print('number of trajectories = ', np.shape(Y)[0])
print('number of frames = ', np.shape(Y)[1])
print('number of TICA dimensions = ',np.shape(Y)[2])

cluster = pyemma.load('cluster.pyemma', model_name='cluster')
dtrajs = cluster.dtrajs

its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=5, errors='bayes')
mplt.plot_implied_timescales(its, outfile="msm_its.png",dt=10, units='ps');
plt.clf()

msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=10, dt_traj='10 ps')
mplt.plot_cktest(msm.cktest(nstates), dt=10, units='ps');
plt.savefig("msm_ck.png")
plt.clf()
msm.save('msm.pyemma', model_name='msm', overwrite=True)
print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))
print('sum of weights = {:f}'.format(msm.pi.sum()))

msm.pcca(nstates)
pcca_samples = msm.sample_by_distributions(msm.metastable_distributions, 10)
pyemma.coordinates.save_trajs(
    reader,
    pcca_samples,
    outfiles=['./data/pcca{}_10samples.pdb'.format(n + 1)
              for n in range(msm.n_metastable)])

bayesian_msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=10, conf=0.95, dt_traj='10 ps')
mplt.plot_cktest(bayesian_msm.cktest(nstates), dt=10, units='ps');
plt.savefig("bayesianmsm_ck.png")
plt.clf()
bayesian_msm.save('bayesian_msm.pyemma', model_name='bayesian_msm', overwrite=True)

hmm = msm.coarse_grain(4)
mplt.plot_markov_model(hmm)
plt.savefig("markov.png")

matplotlib.rcParams.update({'font.size': 14})
dt = 0.1
plt.figure(figsize=(8,5))
ax1=plt.subplot(311)
x = dt*np.arange(tica_traj[0].shape[0])
plt.plot(x, tica_traj[0][:,0]); plt.ylabel('IC 1'); plt.xticks([]); plt.yticks(np.arange(-8, 4, 2))
ax1=plt.subplot(312)
plt.plot(x, tica_traj[0][:,1]); plt.ylabel('IC 2'); plt.xticks([]);  plt.yticks(np.arange(-6, 4, 2))
ax1=plt.subplot(313)
plt.plot(x, tica_traj[0][:,2]); plt.xlabel('time / ns'); plt.ylabel('IC 3'); plt.yticks(np.arange(-4, 6, 2));
plt.savefig("tica_traj.png")
plt.clf()

mplt.plot_free_energy(*tica_traj.T[0:2])
cc_x = cluster.clustercenters[:,0]
cc_y = cluster.clustercenters[:,1]
plt.plot(cc_x,cc_y, linewidth=0, marker='o', markersize=5, color='black')
plt.savefig("cluster.png")
plt.clf()

# xall = np.vstack(tica_traj)[:,0]
# yall = np.vstack(tica_traj)[:,1]
weights=np.concatenate(msm.trajectory_weights())
mplt.plot_free_energy(*tica_traj.T[0:2])
plt.savefig("msm_freeenergy.png")
plt.clf()

proj_ev_all = [np.hstack([msm.eigenvectors_right()[:,i][dtraj] for dtraj in msm.discrete_trajectories_full])
               for i in range(1, 10)]
fig, axes = plt.subplots(1, 3, figsize=(16,4))
for i, ax in enumerate(axes):
    plot_sampled_function(xall, yall, proj_ev_all[i], ax=ax, cbar=False, cmap=plt.cm.Blues)
plt.savefig("eigenvector.png")
plt.clf()

mfpt = np.zeros((nstates, nstates))
for i in range(nstates):
    for j in range(nstates):
        mfpt[i, j] = msm.mfpt(
            msm.metastable_sets[i],
            msm.metastable_sets[j])

inverse_mfpt = np.zeros_like(mfpt)
nz = mfpt.nonzero()
inverse_mfpt[nz] = 1.0 / mfpt[nz]

mplt.plot_network(
    inverse_mfpt,
    pos=np.asarray([[0, 0], [4, 0], [2, 4], [6, 4]]),
    arrow_label_format='%.1f fs',
    arrow_labels=mfpt,
    arrow_scale=3.0,
    state_labels=range(1, nstates + 1),
    size=12);
plt.savefig("network.png")
plt.clf()

A = msm.metastable_sets[0]
B = msm.metastable_sets[1]
flux = pyemma.msm.tpt(msm, A, B)
cg, cgflux = flux.coarse_grain(msm.metastable_sets)
highest_membership = msm.metastable_distributions.argmax(1)
coarse_state_centers = cluster.clustercenters[msm.active_set[highest_membership]]
mplt.plot_flux(
    cgflux,
    coarse_state_centers,
    cgflux.stationary_distribution,
    show_committor=False,
    figpadding=0.2,
    state_labels=['A', '', '', 'B'],
    arrow_label_format='%2.e / ps')
plt.savefig("flux.png")
plt.clf()