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
import glob

topfile=sys.argv[1]
trajfile_path = sys.argv[2]
# Generate a list of files with the "*.nc" format in the directory
trajfile = glob.glob(trajfile_path + '/*_prod.combine.nc')
print(trajfile)
prefix=sys.argv[3]
nstates = 8
nstates1 = 5
nsave_trajs = 100
ncktest = 5
nclusters = 250
nlag = 10

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

tica_ref = pyemma.load('new_tica.pyemma', model_name="tica")

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

data = pyemma.coordinates.load(trajfile, features=feat)

tica = tica_ref.transform(data)
tica_traj = tica # get tica coordinates
tica_concatenated = np.concatenate(tica_traj)
cluster = pyemma.coordinates.cluster_kmeans(tica, k=nclusters, max_iter=50)
dtrajs = cluster.dtrajs
dtrajs_concatenated = np.concatenate(cluster.dtrajs)
#implied timescales (ITS)
its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=5, errors='bayes')
mplt.plot_implied_timescales(its, outfile=prefix+"msm_its.png",dt=10, units='ps');
plt.clf()

msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=nlag, dt_traj='10 ps')

#Chapman-Kolmogorov test
mplt.plot_cktest(msm.cktest(ncktest), dt=10, units='ps');
plt.savefig(prefix+"msm_ck.png")
plt.clf()
msm.save('msm.pyemma', model_name='msm', overwrite=True)

print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))
print('sum of weights = {:f}'.format(msm.pi.sum()))

bayesian_msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=nlag, conf=0.95, dt_traj='10 ps')
mplt.plot_cktest(bayesian_msm.cktest(ncktest), dt=10, units='ps');
plt.savefig(prefix+"bayesianmsm_ck.png")
plt.clf()
bayesian_msm.save('bayesian_msm.pyemma', model_name='bayesian_msm', overwrite=True)

fig, ax = plt.subplots(figsize=(4, 4))
pyemma.plots.plot_density(
    *tica_concatenated[:, :2].T, ax=ax, cbar=False, alpha=0.3)
ax.scatter(*cluster.clustercenters[:, :2].T, s=5, c='C1')
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')
fig.tight_layout()
plt.savefig(prefix+"tica_cluster.png")
plt.clf()

#TICA projection

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

#stationary distribution and the free energy computed over the first two TICA coordinates

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
mplt.plot_contour(
    *tica_concatenated[:, :2].T,
    msm.pi[dtrajs_concatenated],
    ax=axes[0],
    mask=True,
    cbar_label='stationary distribution')
mplt.plot_free_energy(
    *tica_concatenated[:, :2].T,
    weights=np.concatenate(msm.trajectory_weights()),
    ax=axes[1],
    legacy=False)
for ax in axes.flat:
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
axes[0].set_title('Stationary distribution', fontweight='bold')
axes[1].set_title('Reweighted free energy surface', fontweight='bold')
fig.tight_layout()
plt.savefig(prefix+"freeenergy.png")

#eigenvectors corresponding to the slowest processes (largest implied timescales)
eigvec = msm.eigenvectors_right()
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
plt.savefig(prefix+"eigenvectorIC12.png")
plt.clf()

fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_concatenated[:, 2:4].T,
        eigvec[dtrajs_concatenated, i + 1],
        ax=ax,
        cmap='PiYG',
        cbar_label='{}. right eigenvector'.format(i + 2),
        mask=True)
    ax.set_xlabel('IC 3')
axes[0].set_ylabel('IC 4')
fig.tight_layout()
plt.savefig(prefix+"eigenvectorIC34.png")
plt.clf()

fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_concatenated[:, 4:6].T,
        eigvec[dtrajs_concatenated, i + 1],
        ax=ax,
        cmap='PiYG',
        cbar_label='{}. right eigenvector'.format(i + 2),
        mask=True)
    ax.set_xlabel('IC 5')
axes[0].set_ylabel('IC 6')
fig.tight_layout()
plt.savefig(prefix+"eigenvectorIC56.png")
plt.clf()

def coarse_states(nstates, msm, prefix):
    prefix=prefix+str(nstates)
    msm.pcca(nstates)
    pcca_samples = msm.sample_by_distributions(msm.metastable_distributions, nsave_trajs)
    reader = pyemma.coordinates.source(trajfile, features=feat)
    pyemma.coordinates.save_trajs(
        reader,
        pcca_samples,
        outfiles=['./data/pcca{}_{}samples.pdb'.format(n + 1, nsave_trajs)
                for n in range(msm.n_metastable)])

    hmm = msm.coarse_grain(nstates)
    mplt.plot_markov_model(hmm)
    plt.savefig(prefix+"markov.png")

    metastable_traj = msm.metastable_assignments[dtrajs_concatenated]

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
            mfpt[i, j] = msm.mfpt(
                msm.metastable_sets[i],
                msm.metastable_sets[j])
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
        arrow_label_format='%2.e / ps')
    plt.savefig(prefix+"flux.png")
    plt.clf()

coarse_states(nstates, msm, prefix)
coarse_states(nstates1, msm, prefix)
prefix=prefix+"bayesian"
coarse_states(nstates, bayesian_msm, prefix)
coarse_states(nstates1, bayesian_msm, prefix)