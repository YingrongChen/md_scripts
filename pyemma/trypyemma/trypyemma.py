import matplotlib.pyplot as plt
import numpy as np
import pyemma
import sys

topfile=('S129_DNEAY_5ni9.pdb')
train_files=('S129_DNEAY_5ni9_prod.combine.nc')
# test_file=('S129_DNEAY_5ni9_trial1_prod.offset_combine.nc')
feat = pyemma.coordinates.featurizer(topfile)

# feat.add_selection(feat.select_Backbone())
# data_Backbone_position = pyemma.coordinates.load(train_files, features=feat)
indices=[869,871,878,879,880,882,986,987,988,990,1000,1001,1046,1083,1089,1090,1091,1093,2507,2611,2612,2615,2633,2688,2946,2947,2948,2950,2960,2961,2962,2964,2975,2976,2977,2979,2985,2986,2987,2989,3006,3007,3008,3010,3021,3022,3023,3025,3038,3039,3040,3050,3052,3053,3054,3056,3066,3067,3068,3070,3081,3082,3083,3085,3096,3097,3098,3100,3103,3104,3105,3107,3124,3125,3126]
feat.add_selection(feat.add_distances(indices, periodic=False))
data = pyemma.coordinates.load(train_files, features=feat)
print('data dimension: traj time steps, features', data.shape)
#T is the number of time steps in the trajectory 
#d is the number of features (coordinates, observables)
# vamp_Ca_position = pyemma.coordinates.vamp(data_Ca_position, lag=10, dim=3).score(
#         test_data=data_Ca_position,
#         score_method='VAMP2')
# print(f'VAMP2-score xyz, lag=10: {vamp_Ca_position}')

tica = pyemma.coordinates.tica(data)
tica_concatenated = np.concatenate(tica.get_output())
print('tica dimension: traj time steps, tica', tica_concatenated.shape)

# cluster in tica space
cluster = pyemma.coordinates.cluster_kmeans(tica, k=100, max_iter=50)
# print('cluster dimension:', cluster.describe())

# pyemma.plots.plot_feature_histograms(
#     tica_concatenated, ['IC {}'.format(i + 1) for i in range(tica.dimension())])
# plt.savefig('tica_feature.png')

its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10], nits=3, errors='bayes')
pyemma.plots.plot_implied_timescales(its, ylog=False);
plt.savefig('msm_its.png')

msm = pyemma.msm.estimate_markov_model(cluster.dtrajs, lag=1)
pyemma.plots.plot_cktest(msm.cktest(2));
plt.savefig('msm_estimate_markov_model.png')

# vamp_Ca_position = pyemma.coordinates.vamp(data_Ca_position, lag=2, dim=3).score(
#         test_data=data_Ca_position,
#         score_method='VAMP2')
# print(f'VAMP2-score xyz, lag=2: {vamp_Ca_position}')

# data_torsions_test = pyemma.coordinates.load(test_file, features=feat)

# feat.active_features = []
# feat.add_distances_ca(periodic=False)
# data_dists_ca = pyemma.coordinates.load(train_files, features=feat)
# # data_dists_ca_test = pyemma.coordinates.load(test_file, features=feat)

# feat.active_features = []
# pairs = feat.pairs(feat.select_Heavy())
# feat.add_contacts(pairs, periodic=False)
# data_contacts = pyemma.coordinates.load(train_files, features=feat)
# # data_contacts_test = pyemma.coordinates.load(test_file, features=feat)

# vamp_dist_ca = pyemma.coordinates.vamp(data_dists_ca, lag=lag, dim=dim)
# vamp_contacts = pyemma.coordinates.vamp(data_contacts, lag=lag, dim=dim)
# print('VAMP2-score distances between all Ca carbon atoms: {:f}'.format(vamp_dist_ca))
# print('VAMP2-score xyz: {:f}'.format(vamp_contacts))

# def plot_for_lag(ax, lag, dim=3):
#     # vamp_torsions = pyemma.coordinates.vamp(data_torsions, lag=lag, dim=dim)
#     vamp_dist_ca = pyemma.coordinates.vamp(data_dists_ca, lag=lag, dim=dim)
#     vamp_contacts = pyemma.coordinates.vamp(data_contacts, lag=lag, dim=dim)
#     print('VAMP2-score backbone torsions: {:f}'.format(vamp_dist_ca))
#     print('VAMP2-score xyz: {:f}'.format(vamp_contacts))
#     # vamps = (vamp_dist_ca, vamp_contacts)
#     # test_data = (data_dists_ca_test, data_contacts_test)
#     # labels = ('CA distances', 'contacts')
#     # for i, (v, test_data) in enumerate(zip(vamps, test_data)):
#     #     s = v.score(test_data=test_data)
#     #     ax.bar(i, s)
#     # ax.set_title('VAMP2 @ lag = {} ps'.format(lag))
#     # ax.set_xticks(range(len(vamps)))
#     # ax.set_xticklabels(labels)
#     # fig.tight_layout()
#     # plt.savefig(str(lag),'.png')

# # fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharey=True)
# plot_for_lag(axes[0], 5)
# plot_for_lag(axes[1], 10)
# plot_for_lag(axes[2], 20)
# plot_for_lag(axes[3], 50)