#export PATH=/dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
#setenv PATH /dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH

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
from pandas import DataFrame

topfile=sys.argv[1]
trajfile=sys.argv[2]
prefix=sys.argv[3]
nstates = 5
nsave_trajs = 20
nclusters = 300
nlag = 20

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

def featurize_data(topfile):
    feat = pyemma.coordinates.featurizer(topfile)
    ser=feat.select("(resSeq 189) and name CA")
    tyr=feat.select("(resSeq 185) and name CA")
    # nottail="(resSeq 4 to 81 and 88 to 175)"
    notser=feat.select("(not resSeq 189) and name CA")
    nottyr=feat.select("(not (resSeq 189 or resSeq 185)) and name CA")
    ahelix = slice(51, 77)
    bhelix = slice(139, 172)
    peptide = slice(181, 194)
    n = len(range(ahelix.start, ahelix.stop))
    m = len(range(peptide.start, peptide.stop))
    matrix = np.zeros((3 * n + 6 * m, 2), dtype=int)
    for j in range(0, 3):
        for i in range(0, n):
            matrix[i + j * n, :] = [bhelix.start + i + j, ahelix.stop - 1 - i]
    for j in range(0, 3):
        for i in range(0, m):
            matrix[i + 3 * n + j * m, :] = [peptide.start + i, ahelix.stop - 1 - i - j * 2]
    for j in range(0, 3):
        for i in range(0, m):
            matrix[i + 3 * n + (3 + j) * m, :] = [peptide.start + i, bhelix.stop - 1 - i - j * 2]

    feat.add_backbone_torsions()
    feat.add_distances(notser, indices2=ser)
    feat.add_distances(nottyr, indices2=tyr)
    feat.add_residue_mindist(residue_pairs=matrix)
    print('number of features ', feat.dimension())
    return feat

feat = featurize_data(topfile)
# reader = pyemma.coordinates.source(trajfile, features=feat)
data = pyemma.coordinates.load(trajfile, features=feat)

if os.path.exists(prefix+'tica.pyemma'):
    print("TICA exists!")
    tica = pyemma.load(prefix+'tica.pyemma', model_name='tica')
    tica_traj = tica.transform(data)
    tica_concatenated = tica_traj
else:
    print("TICA does not exist.")
    tica = pyemma.coordinates.tica(data, lag=nlag)
    tica.save(prefix+'tica.pyemma', model_name='tica', overwrite=True)
    tica_traj = tica.get_output() # get tica coordinates
    tica_concatenated = np.concatenate(tica_traj)

    matplotlib.rcParams.update({'font.size': 14})
    dt = 0.1
    plt.figure(figsize=(10,5))
    ax1=plt.subplot(511)
    x = dt*np.arange(tica_traj[0].shape[0])
    plt.plot(x, tica_traj[0][:,0]); plt.ylabel('IC 1'); plt.xticks([]); plt.yticks(np.arange(-8, 4, 2))
    ax1=plt.subplot(512)
    plt.plot(x, tica_traj[0][:,1]); plt.ylabel('IC 2'); plt.xticks([]);  plt.yticks(np.arange(-6, 4, 2))
    ax1=plt.subplot(513)
    plt.plot(x, tica_traj[0][:,2]); plt.ylabel('IC 3'); plt.yticks(np.arange(-4, 6, 2));
    plt.savefig(prefix+"tica_traj.png")
    ax1=plt.subplot(514)
    plt.plot(x, tica_traj[0][:,3]); plt.ylabel('IC 4'); plt.xticks([]);  plt.yticks(np.arange(-6, 4, 2))
    ax1=plt.subplot(515)
    plt.plot(x, tica_traj[0][:,4]); plt.xlabel('time / ns'); plt.ylabel('IC 5'); plt.yticks(np.arange(-4, 6, 2));
    plt.savefig(prefix+"tica_traj.png")
    plt.clf()

print('shape of TICA = ',np.shape(tica_traj)) #(400000, 321)

if os.path.exists(prefix+'cluster.pyemma'):
    print("Cluster exists!")
    cluster = pyemma.load(prefix+'cluster.pyemma', model_name='cluster')
else:
    print("Cluster does not exist.")
    cluster = coor.cluster_kmeans(tica_traj, k=nclusters, max_iter=50)
    cluster.save(prefix+'cluster.pyemma', model_name='cluster', overwrite=True)
    fig, ax = plt.subplots(figsize=(4, 4))

    pyemma.plots.plot_density(
        *tica_concatenated[:, :2].T, ax=ax, cbar=False, alpha=0.3)
    ax.scatter(*cluster.clustercenters[:, :2].T, s=5, c='C1')
    ax.set_xlabel('IC 1')
    ax.set_ylabel('IC 2')
    fig.tight_layout()
    plt.savefig(prefix+"tica_cluster.png")
    plt.clf()

dtrajs_concatenated = np.concatenate(cluster.dtrajs)
print('shape of tica_concatenated = ',np.shape(tica_concatenated)) #(400000, 321)

# #implied timescales (ITS)
# its = pyemma.msm.its(cluster.dtrajs, lags=[1, 5, 10, 15, 20, 25, 30], nits=5, errors='bayes')
# mplt.plot_implied_timescales(its, outfile=prefix+"msm_its.png",dt=10, units='ps');
# plt.clf()

# its_hmm = pyemma.msm.timescales_hmsm(cluster.dtrajs, nstates, lags=[1, 5, 10, 15, 20, 25, 30], nits=5, errors='bayes')
# mplt.plot_implied_timescales(its_hmm, outfile="hmm_its.png",dt=10, units='ps');
# plt.clf()

# msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=nlag, dt_traj='10 ps')
# msm.save(prefix+'msm.pyemma', model_name='msm', overwrite=True)

hmm = pyemma.msm.estimate_hidden_markov_model(cluster.dtrajs, lag=nlag, nstates=nstates, dt_traj='10 ps')
hmm.save(prefix+'hmm.pyemma', model_name='hmm', overwrite=True)
#Chapman-Kolmogorov test
mplt.plot_cktest(hmm.cktest(nstates), dt=10, units='ps');
plt.savefig(prefix+"hmm_ck.png")
plt.clf()

print("transition_matrix, stationary_distribution, observable_set, observation_probabilities, lifetimes, timescales, mfpt")
print(hmm.transition_matrix)
print(hmm.stationary_distribution)
print(hmm.observable_set)
print(hmm.observation_probabilities) 
print(hmm.lifetimes) 
print(hmm.timescales())  

#stationary distribution and the free energy computed over the first two TICA coordinates
fig, ax, misc = mplt.plot_free_energy(
    *tica_concatenated[:, :2].T,
    weights=np.concatenate(hmm.trajectory_weights()),
    legacy=False)
fig.tight_layout()
plt.savefig(prefix+"freeenergy.png")

#eigenvectors corresponding to the slowest processes (largest implied timescales)

eigvec = hmm.eigenvectors_right()
print('The first eigenvector is one: {} (min={}, max={})'.format(
    np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))
fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_concatenated[:, :2].T,
        eigvec[dtrajs_concatenated, i + 1],
        ax=ax,
        cmap='PiYG',
        cbar_label='{}. right eigenvector'.format(i + 2),
        mask=True)
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
fig.tight_layout()
plt.savefig(prefix+"eigenvector.png")
plt.clf()

pcca_samples = hmm.sample_by_distributions(hmm.metastable_distributions, nsave_trajs)
reader = pyemma.coordinates.source(trajfile, features=feat)
pyemma.coordinates.save_trajs(
    reader,
    pcca_samples,
    outfiles=['./data/pcca{}_{}samples.pdb'.format(n + 1, nsave_trajs)
            for n in range(hmm.n_metastable)])

metastable_traj = hmm.metastable_assignments[dtrajs_concatenated]

fig, ax = plt.subplots(figsize=(5, 4))
_, _, misc = pyemma.plots.plot_state_map(
    *tica_concatenated[:, :2].T, metastable_traj, ax=ax)
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')
misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
                            for i in range(nstates)])
fig.tight_layout()
plt.savefig(prefix+"memberships.png")
plt.clf()

mfpt = np.zeros((nstates, nstates)) #mean first passage times (MFPTs)
for i in range(nstates):
    for j in range(nstates):
        mfpt[i, j] = hmm.mfpt(
            hmm.metastable_sets[i],
            hmm.metastable_sets[j])
df = DataFrame(np.round(mfpt, decimals=2), index=range(1, nstates + 1), columns=range(1, nstates + 1))
df.to_csv(prefix+'MFPT.csv')  

inverse_mfpt = np.zeros_like(mfpt)
nz = mfpt.nonzero()
inverse_mfpt[nz] = 1.0 / mfpt[nz]

mplt.plot_network(
    inverse_mfpt,
    arrow_label_format='%.1f fs',
    arrow_labels=mfpt,
    arrow_scale=3.0,
    state_labels=range(1, nstates + 1),
    size=12);
plt.savefig(prefix+"network.png")
plt.clf()

A = hmm.metastable_sets[0]
B = hmm.metastable_sets[1]
flux = pyemma.hmm.tpt(hmm, A, B)
cg, cgflux = flux.coarse_grain(hmm.metastable_sets)
highest_membership = hmm.metastable_distributions.argmax(1)
coarse_state_centers = cluster.clustercenters[hmm.active_set[highest_membership]]
mplt.plot_flux(
    cgflux,
    coarse_state_centers,
    cgflux.stationary_distribution,
    show_committor=False,
    figpadding=0.2,
    arrow_label_format='%2.e / ps')
plt.savefig(prefix+"flux.png")
plt.clf()