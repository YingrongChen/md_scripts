#construct a shared backbone tica model for all the trajectories and map individual trajectories to the shared tica space
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
ref_trajfile = '../project/bb_nontcr_4x5w_5ni9.nc'
trajfile = sys.argv[2]
prefix = sys.argv[3]
traj_size = int(sys.argv[4])
nclusters = 250
ncktest = 5
nlag = 10
nstates = 3
nstates1 = 5
nsave_trajs = 100
restart=2

def featurize_data(topfile):
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
    return feat

with open("succes_concate.txt", 'r') as file:
    segment_names = [line.strip().split()[0] for line in file.readlines()]
segment_pointer = restart*5

feat = featurize_data(topfile)
# reader = pyemma.coordinates.source(trajfile, features=feat)
ref_data = pyemma.coordinates.load(ref_trajfile, features=feat)
ref_tica = pyemma.load('ref_tica.pyemma', model_name='tica')
ref_tica_concatenated = ref_tica.transform(ref_data)

data = pyemma.coordinates.load(trajfile, features=feat)
combine_tica_concatenated = ref_tica.transform(data)

print('shape of reference TICA = ',np.shape(ref_tica_concatenated)) #(400000, 321)
print('shape of combine TICA = ',np.shape(combine_tica_concatenated)) #(400000, 321)
print(f"Number of names in the list { len(segment_names) // 5 } ")
print(f"should match the number of trajectory { np.shape(combine_tica_concatenated)[0] // (traj_size * 5) }")


for i in range(restart*traj_size, np.shape(combine_tica_concatenated)[0], traj_size * 5):
    segment_name = segment_names[segment_pointer]
    segment_pointer = segment_pointer + 5
    print(segment_name)
    # Create a directory for the segment if it does not exist
    segment_dir = os.path.join(os.getcwd(), segment_name)
    os.makedirs(segment_dir, exist_ok=True)    
    
    tica_traj = []
    for j in range(5):
        segment = combine_tica_concatenated[i + j * traj_size: i + (j + 1) * traj_size]
        tica_traj.append(segment)
    tica_concatenated = np.concatenate(tica_traj)
    cluster = coor.cluster_kmeans(tica_traj, k=nclusters, max_iter=50)
    dtrajs_concatenated = np.concatenate(cluster.dtrajs)
    cluster.save(os.path.join(segment_dir, 'cluster.pyemma'), model_name='cluster', overwrite=True)

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
    plt.savefig(os.path.join(segment_dir, "tica_cluster.png"))
    plt.clf()
    
    its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=5, errors='bayes')
    mplt.plot_implied_timescales(its, outfile=os.path.join(segment_dir,"msm_its.png"),dt=10, units='ps');
    plt.clf()

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
    plt.savefig(os.path.join(segment_dir,"tica_traj.png"))
    plt.clf()