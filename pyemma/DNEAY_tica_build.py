#export PATH=/dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH
#setenv PATH /dors/meilerlab/home/brownbp1/miniconda3/envs/pyemma/bin:$PATH

import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
import warnings

# Set the paths to the top and trajectory files
topfile = sys.argv[1]
trajfile = sys.argv[2]
prefix = "S_"
traj_size = 100000
nclusters = 250
ncktest = 5
nlag = 10
nstates = 3
nstates1 = 5
nsave_trajs = 100
restart=1

segment_names = ["DNEAY_5ni9", "S129_DNEAY_5ni9"]

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
ref_data = pyemma.coordinates.load(trajfile, features=feat)

if os.path.exists(prefix+'ref_tica.pyemma'):
    print("TICA exists!")
    ref_tica = pyemma.load(prefix+'ref_tica.pyemma', model_name='tica')
    ref_tica_concatenated = ref_tica.transform(ref_data)
else:
    print("TICA does not exist.")
    ref_tica = pyemma.coordinates.tica(ref_data, lag=nlag)
    ref_tica.save(prefix+'ref_tica.pyemma', model_name='tica', overwrite=True)
    ref_tica_traj = ref_tica.get_output() # get ref_tica coordinates
    ref_tica_concatenated = np.concatenate(ref_tica_traj)

print('shape of reference TICA = ',np.shape(ref_tica_concatenated)) #(400000, 321)
print(f"Number of names in the list { len(segment_names) // 5 } ")
print(f"should match the number of trajectory { np.shape(ref_tica_concatenated)[0] // (traj_size * 5) }")
segment_pointer=restart
for i in range(restart*traj_size, np.shape(ref_tica_concatenated)[0], traj_size * 5):
    segment_name = segment_names[segment_pointer]
    segment_pointer = segment_pointer + 1
    print(segment_name)

    segment_dir = os.path.join(os.getcwd(), segment_name)
    os.makedirs(segment_dir, exist_ok=True)    

    tica_traj = []
    for j in range(5):
        segment = ref_tica_concatenated[i + j * traj_size: i + (j + 1) * traj_size]
        tica_traj.append(segment)
    tica_concatenated = np.concatenate(tica_traj)


# if os.path.exists(os.path.join(segment_dir, prefix+'cluster.pyemma')):
#     print("Cluster exists!")
#     cluster = pyemma.load(os.path.join(segment_dir, prefix+'cluster.pyemma'), model_name='cluster')
#     ref_tica_concatenated = ref_tica.transform(ref_data)
# else:
#     print("Cluster does not exist.")

    cluster = coor.cluster_kmeans(tica_traj, k=nclusters, max_iter=50)
    cluster.save(os.path.join(segment_dir, prefix+'cluster.pyemma'), model_name='cluster', overwrite=True)

    dtrajs_concatenated = np.concatenate(cluster.dtrajs)
    
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
    plt.savefig(os.path.join(segment_dir, prefix+"tica_cluster.png"))
    plt.clf()
    
    #implied timescales (ITS)
    # its = pyemma.msm.its(cluster.dtrajs, lags=[1,5,10,15,20,25,30], nits=5, errors='bayes')
    # mplt.plot_implied_timescales(its, outfile=os.path.join(segment_dir,prefix+"msm_its.png"),dt=10, units='ps');
    # plt.clf()

    # its_hmm = pyemma.msm.timescales_hmsm(cluster.dtrajs, nstates, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=5, errors='bayes')
    # mplt.plot_implied_timescales(its_hmm, outfile=os.path.join(segment_dir,"hmm_its.png"),dt=10, units='ps');
    # plt.clf()

    # # msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=nlag, dt_traj='10 ps')
    # hmm = pyemma.msm.estimate_hidden_markov_model(cluster.dtrajs, nstates, lag=nlag, dt_traj='10 ps')
    # #Chapman-Kolmogorov test
    # mplt.plot_cktest(hmm.cktest(ncktest), dt=10, units='ps');
    # plt.savefig(os.path.join(segment_dir, prefix+"hmm_ck.png"))
    # plt.clf()
    # hmm.save(os.path.join(segment_dir, prefix+'hmm.pyemma'), model_name='hmm', overwrite=True)

    dt = 0.1
    plt.figure(figsize=(10,5))
    ax1=plt.subplot(511)
    x = dt*np.arange(tica_concatenated.shape[0])
    plt.plot(x, tica_concatenated[:,0]); plt.ylabel('IC 1'); plt.xticks([]); plt.yticks(np.arange(-8, 4, 2))
    ax1=plt.subplot(512)
    plt.plot(x, tica_concatenated[:,1]); plt.ylabel('IC 2'); plt.xticks([]);  plt.yticks(np.arange(-6, 4, 2))
    ax1=plt.subplot(513)
    plt.plot(x, tica_concatenated[:,2]); plt.ylabel('IC 3'); plt.xticks([]);  plt.yticks(np.arange(-4, 4, 2))
    ax1=plt.subplot(514)
    plt.plot(x, tica_concatenated[:,3]); plt.ylabel('IC 4'); plt.xticks([]);  plt.yticks(np.arange(-2, 4, 2))
    ax1=plt.subplot(515)
    plt.plot(x, tica_concatenated[:,4]); plt.xlabel('time / ns'); plt.ylabel('IC 5'); plt.yticks(np.arange(0, 6, 2))
    plt.savefig(os.path.join(segment_dir,prefix+"tica_traj.png"))
    plt.clf()

    # # #stationary distribution and the free energy computed over the first two TICA coordinates
    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    # mplt.plot_free_energy(
    #     *tica_concatenated[:, [0, 1]].T,
    #     weights=np.concatenate(msm.trajectory_weights()),
    #     ax=ax1,
    #     legacy=False)
    # mplt.plot_free_energy(
    #     *tica_concatenated[:, [2, 3]].T,
    #     weights=np.concatenate(msm.trajectory_weights()),
    #     ax=ax2,
    #     legacy=False)
    # mplt.plot_free_energy(
    #     *tica_concatenated[:, [4, 5]].T,
    #     weights=np.concatenate(msm.trajectory_weights()),
    #     ax=ax3,
    #     legacy=False)
    # fig.tight_layout()
    # plt.savefig(os.path.join(segment_dir,"freeenergy.png"))

    # #eigenvectors corresponding to the slowest processes (largest implied timescales)
    # # eigvec = msm.eigenvectors_right()
    # # print('The first eigenvector is one: {} (min={}, max={})'.format(
    # #     np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))
    # # fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
    # # for j, ax in enumerate(axes.flat):
    # #     pyemma.plots.plot_contour(
    # #         *tica_concatenated[:, :2].T,
    # #         eigvec[dtrajs_concatenated, j + 1],
    # #         ax=ax,
    # #         cmap='PiYG',
    # #         cbar_label='{}. right eigenvector'.format(j + 2),
    # #         mask=True)
    # #     ax.set_xlabel('IC 1')
    # # axes[0].set_ylabel('IC 2')
    # # fig.tight_layout()
    # # plt.savefig(os.path.join(segment_dir,"eigenvectorIC12.png"))
    # # plt.clf()

    # msm.pcca(nstates)
    # pcca_samples = msm.sample_by_distributions(msm.metastable_distributions, nsave_trajs)
    # pyemma.coordinates.save_trajs(
    #     reader,
    #     pcca_samples,
    #     outfiles = [os.path.join(segment_dir, 'pcca{}_{}samples.pdb'.format(n + 1, nsave_trajs)) for n in range(msm.n_metastable)])

    # metastable_traj = msm.metastable_assignments[dtrajs_concatenated]

    # fig, ax = plt.subplots(figsize=(5, 4))
    # _, _, misc = pyemma.plots.plot_state_map(
    #     *tica_concatenated[:, :2].T, metastable_traj, ax=ax)
    # ax.set_xlabel('IC 1')
    # ax.set_ylabel('IC 2')
    # misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
    #                             for i in range(nstates)])
    # fig.tight_layout()
    # plt.savefig(os.path.join(segment_dir, "memberships.png"))
    # plt.clf()

    # mfpt = np.zeros((nstates, nstates)) #mean first passage times (MFPTs)
    # for i in range(nstates):
    #     for j in range(nstates):
    #         mfpt[i, j] = msm.mfpt(
    #             msm.metastable_sets[i],
    #             msm.metastable_sets[j])
    # df = DataFrame(np.round(mfpt, decimals=2), index=range(1, nstates + 1), columns=range(1, nstates + 1))
    # df.to_csv(os.path.join(segment_dir, 'MFPT.csv'))

    # inverse_mfpt = np.zeros_like(mfpt)
    # nz = mfpt.nonzero()
    # inverse_mfpt[nz] = 1.0 / mfpt[nz]

    # A = msm.metastable_sets[0]
    # B = msm.metastable_sets[1]
    # flux = pyemma.msm.tpt(msm, A, B)
    # cg, cgflux = flux.coarse_grain(msm.metastable_sets)
    # highest_membership = msm.metastable_distributions.argmax(1)
    # coarse_state_centers = cluster.clustercenters[msm.active_set[highest_membership]]
    # mplt.plot_flux(
    #     cgflux,
    #     coarse_state_centers,
    #     cgflux.stationary_distribution,
    #     show_committor=False,
    #     figpadding=0.2,
    #     arrow_label_format='%2.e / ps')
    # plt.savefig(os.path.join(segment_dir, "flux.png"))
    # plt.clf()
