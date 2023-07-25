#export PATH=/dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
#setenv PATH /dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
import glob
import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import warnings
from pandas import DataFrame

topfile = sys.argv[1]
# trajfile = sys.argv[2]
# prefix = sys.argv[3]
# traj_size = int(sys.argv[4])
# prefix=os.path.basename(os.getcwd()) + "_"
prefix=topfile.replace(".pdb", "_")
trajfile=glob.glob(prefix + "trial?_prod.combine.nc")
print(trajfile)
nclusters = 200
ncktest = 5
nlag = 1
nstates = 5
# nstates1 = 2
nsave_trajs = 100
nstride=1

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
    ser=feat.select("(resSeq 188) and name CA")
    tyr=feat.select("(resSeq 184) and name CA")
    tyr2=feat.select("(resSeq 184) and name CZ")
    # nottail="(resSeq 4 to 81 and 88 to 175)"
    notser=feat.select("(not resSeq 188) and name CA")
    nottyr=feat.select("(not (resSeq 188 or resSeq 184)) and name CA")
    # nottyr2=feat.select("(not (resSeq 188 or resSeq 184 or resSeq 192)) and name CA")

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
    feat.add_distances(nottyr, indices2=tyr2)
    feat.add_distances(matrix)
    print('number of features ', feat.dimension())
    return feat

feat = featurize_data(topfile)
data = pyemma.coordinates.load(trajfile, features=feat)
reader = pyemma.coordinates.source(trajfile, features=feat)

# ref_tica = pyemma.load('../ref_tica.pyemma', model_name='tica')
# tica_traj = ref_tica.transform(data)

if os.path.exists(prefix+'tica.pyemma'):
    print("TICA exists!")
    tica = pyemma.load(prefix+'tica.pyemma', model_name='tica')
    tica_traj = tica.transform(data)
    tica_concatenated = np.concatenate(tica_traj)
else:
    print("TICA does not exist.")
    tica = pyemma.coordinates.tica(data)
    tica_traj = tica.get_output() # get ref_tica coordinates
    tica_concatenated = np.concatenate(tica_traj)
    tica.save(prefix+'tica.pyemma', model_name='tica', overwrite=True)

# if os.path.exists(prefix+'msm.pyemma'):
#     print("MSM exists!")
#     cluster = pyemma.load(prefix+'cluster.pyemma', model_name='cluster')
#     msm = pyemma.load(prefix+'msm.pyemma', model_name='msm')
# else:
cluster = coor.cluster_kmeans(tica_traj, k=nclusters, max_iter=50)
cluster.save(prefix+'cluster.pyemma', model_name='cluster', overwrite=True)

its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=5, errors='bayes')
mplt.plot_implied_timescales(its, outfile=prefix+"msm_its.png",dt=10, units='ps')
plt.clf()
msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=nlag, dt_traj='10 ps')
msm.save(prefix+'msm.pyemma', model_name='msm', overwrite=True)

mplt.plot_cktest(msm.cktest(ncktest), dt=10, units='ps');
plt.savefig(prefix+"msm_ck.png")
plt.clf()

os.makedirs('data', exist_ok=True)
msm.pcca(nstates)
pcca_samples = msm.sample_by_distributions(msm.metastable_distributions, nsave_trajs)
pyemma.coordinates.save_trajs(
    reader,
    pcca_samples,
    outfiles = ['./data/{}{}samples.pdb'.format(prefix, n + 1) for n in range(msm.n_metastable)])

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

mplt.plot_density(*ref_tica_concatenated[:, [0, 1]].T, ax=ax1, cbar=False, alpha=0.3)
ax1.scatter(*cluster.clustercenters[:, [0, 1]].T, s=5, c='C1', label='Cluster centers')
ax1.set_xlabel('IC 1')
ax1.set_ylabel('IC 2')
ax1.legend()

mplt.plot_density(*ref_tica_concatenated[:, [2, 3]].T, ax=ax2, cbar=False, alpha=0.3)
ax2.scatter(*cluster.clustercenters[:, [2, 3]].T, s=5, c='C1', label='Cluster centers')
ax2.set_xlabel('IC 3')
ax2.set_ylabel('IC 4')
ax2.legend()

mplt.plot_density(*ref_tica_concatenated[:, [4, 5]].T, ax=ax3, cbar=False, alpha=0.3)
ax3.scatter(*cluster.clustercenters[:, [4, 5]].T, s=5, c='C1', label='Cluster centers')
ax3.set_xlabel('IC 5')
ax3.set_ylabel('IC 6')
ax3.legend()

fig.tight_layout()
plt.savefig(prefix+"tica_cluster.png")
plt.clf()

print("tica_traj shape", np.shape(tica_traj))
print("tica_concatenated shape", np.shape(tica_concatenated))
dtrajs_concatenated = np.concatenate(cluster.dtrajs)
print("cluster.dtrajs shape", np.shape(cluster.dtrajs))
print("dtrajs_concatenated shape", np.shape(dtrajs_concatenated))
print('fraction of states used = {:f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:f}'.format(msm.active_count_fraction))
print(msm.active_set)
print(msm.stationary_distribution)
print('sum of weights = {:f}'.format(msm.pi.sum()))

# for i, idx in enumerate(dtrajs_concatenated):
#     if idx not in msm.active_set:
#         print(f"{i}: {idx}")

filtered_indices = [i for i, idx in enumerate(dtrajs_concatenated) if idx in msm.active_set]
dtrajs_filtered = dtrajs_concatenated[filtered_indices]
tica_filtered = tica_concatenated[filtered_indices]

# dtrajs_filtered = dtrajs_concatenated
print("dtrajs_filtered shape", np.shape(dtrajs_filtered))

#stationary distribution and the free energy computed over the first two TICA coordinates
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
mplt.plot_contour(
    *tica_concatenated[:, :2].T,
    msm.pi[dtrajs_filtered],
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
print("eigvec shape", np.shape(eigvec))
print('The first eigenvector is one: {} (min={}, max={})'.format(
    np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))
print('second eigenvector is: min={}, max={}'.format(
    eigvec[:, 1].min(), eigvec[:, 1].max()))
fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
for j, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_concatenated[:, :2].T,
        eigvec[dtrajs_filtered, j + 1],
        ax=ax,
        cmap='PiYG',
        cbar_label='{}. right eigenvector'.format(j + 2),
        mask=True)
    ax.set_xlabel('IC 1')
axes[0].set_ylabel('IC 2')
fig.tight_layout()
plt.savefig(prefix+"eigenvectorIC12.png")
plt.clf()

metastable_traj = msm.metastable_assignments[dtrajs_filtered]
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
