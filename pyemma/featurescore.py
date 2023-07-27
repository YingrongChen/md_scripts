#try different features and see which one is the best

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

scores = []
dimensions = []
distance_feat = pyemma.coordinates.featurizer(topfile)
atoms_heavy = distance_feat.select_Heavy()
ser = distance_feat.select("(resSeq 189 and name CA)")
tyr = distance_feat.select("(resSeq 185 and name CA)")
notser=distance_feat.select("(not resSeq 189) and name CA")
nottyr=distance_feat.select("(not (resSeq 189 or resSeq 185)) and name CA")
distance_feat.add_distances(notser, indices2=ser)
distance_feat.add_distances(nottyr, indices2=tyr)
distance_data = pyemma.coordinates.load(trajfile, features=distance_feat)
score_distance_data = pyemma.coordinates.vamp(distance_data[:-1], dim=2).score(
        test_data=distance_data[:-1],
        score_method='VAMP2')
print(distance_feat.describe())
print('number of single-point-to-others distance features ', distance_feat.dimension())
print('VAMP2-score single-point-to-others distance: {:f}'.format(score_distance_data))

backbone_feat = pyemma.coordinates.featurizer(topfile)
backbone_feat.add_backbone_torsions()
torsions_data = pyemma.coordinates.load(trajfile, features=backbone_feat)
score_torsions_data = pyemma.coordinates.vamp(torsions_data[:-1], dim=2).score(
        test_data=torsions_data[:-1],
        score_method='VAMP2')
print(backbone_feat.describe())
print('number of backbone features ', backbone_feat.dimension())
print('VAMP2-score backbone torsions: {:f}'.format(score_torsions_data))

sidechain_feat = pyemma.coordinates.featurizer(topfile)
sidechain_feat.add_sidechain_torsions(which='chi1')
sidechain_feat.add_sidechain_torsions(which='chi3')
sidechain_data = pyemma.coordinates.load(trajfile, features=sidechain_feat)
score_sidechain_data = pyemma.coordinates.vamp(sidechain_data[:-1], dim=2).score(
        test_data=sidechain_data[:-1],
        score_method='VAMP2')
print(sidechain_feat.describe())
print('number of sidechain features ', sidechain_feat.dimension())
print('VAMP2-score sidechain torsions: {:f}'.format(score_sidechain_data))

pairdistance_feat = pyemma.coordinates.featurizer(topfile)
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
pairdistance_feat.add_residue_mindist(residue_pairs=matrix)
pairdistance_data = pyemma.coordinates.load(trajfile, features=pairdistance_feat)
score_pairdistance_data = pyemma.coordinates.vamp(pairdistance_data[:-1], dim=2).score(
        test_data=pairdistance_data[:-1],
        score_method='VAMP2')
print(pairdistance_feat.describe())
print('number of bhelix-ahelix-peptide pair distance features ', pairdistance_feat.dimension())
print('VAMP2-score bhelix-ahelix-peptide pair distance: {:f}'.format(score_pairdistance_data))

# def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
#     nval = int(len(data) * validation_fraction)
#     scores = np.zeros(number_of_splits)
#     for n in range(number_of_splits):
#         ival = np.random.choice(len(data), size=nval, replace=False)
#         vamp = pyemma.coordinates.vamp(
#             [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
#         scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
#     return scores

# dim = 10
# fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
# for ax, lag in zip(axes.flat, [5, 10, 20]):
#     torsions_scores = score_cv(torsions_data, lag=lag, dim=dim)
#     scores = [torsions_scores.mean()]
#     errors = [torsions_scores.std()]
#     distance_scores = score_cv(distance_data, lag=lag, dim=dim)
#     scores += [distance_scores.mean()]
#     errors += [distance_scores.std()]
#     pairdistance_scores = score_cv(pairdistance_data, lag=lag, dim=dim)
#     scores += [pairdistance_scores.mean()]
#     errors += [pairdistance_scores.std()]
#     sidechain_scores = score_cv(sidechain_data, lag=lag, dim=dim)
#     scores += [sidechain_scores.mean()]
#     errors += [sidechain_scores.std()]
#     ax.bar(labels, scores, yerr=errors, color=['C0', 'C1', 'C2', 'C3'])
#     ax.set_title(r'lag time $\tau$={:.1f}ns'.format(lag * 0.1))
#     if lag == 5:
#         # save for later
#         vamp_bars_plot = dict(
#             labels=labels, scores=scores, errors=errors, dim=dim, lag=lag)
# axes[0].set_ylabel('VAMP2 score')
# fig.tight_layout()
# plt.savefig(prefix+"all_vamp2.png")
# plt.clf()

# lags = [1, 2, 5, 10, 20]
# dims = [i + 1 for i in range(10)]

# fig, ax = plt.subplots()
# for i, lag in enumerate(lags):
#     scores_ = np.array([score_cv(torsions_data, dim, lag)
#                         for dim in dims])
#     scores = np.mean(scores_, axis=1)
#     errors = np.std(scores_, axis=1, ddof=1)
#     color = 'C{}'.format(i)
#     ax.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
#     ax.plot(dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
# ax.legend()
# ax.set_xlabel('number of dimensions')
# ax.set_ylabel('VAMP2 score')
# fig.tight_layout()
# plt.savefig(prefix+"torsions_vamp2.png")
# plt.clf()

# fig, ax = plt.subplots()
# for i, lag in enumerate(lags):
#     scores_ = np.array([score_cv(distance_data, dim, lag)
#                         for dim in dims])
#     scores = np.mean(scores_, axis=1)
#     errors = np.std(scores_, axis=1, ddof=1)
#     color = 'C{}'.format(i)
#     ax.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
#     ax.plot(dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
# ax.legend()
# ax.set_xlabel('number of dimensions')
# ax.set_ylabel('VAMP2 score')
# fig.tight_layout()
# plt.savefig(prefix+"distance_vamp2.png")
# plt.clf()

# fig, ax = plt.subplots()
# for i, lag in enumerate(lags):
#     scores_ = np.array([score_cv(pairdistance_data, dim, lag)
#                         for dim in dims])
#     scores = np.mean(scores_, axis=1)
#     errors = np.std(scores_, axis=1, ddof=1)
#     color = 'C{}'.format(i)
#     ax.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
#     ax.plot(dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
# ax.legend()
# ax.set_xlabel('number of dimensions')
# ax.set_ylabel('VAMP2 score')
# fig.tight_layout()
# plt.savefig(prefix+"pairdistance_vamp2.png")
# plt.clf()

# fig, ax = plt.subplots()
# for i, lag in enumerate(lags):
#     scores_ = np.array([score_cv(sidechain_data, dim, lag)
#                         for dim in dims])
#     scores = np.mean(scores_, axis=1)
#     errors = np.std(scores_, axis=1, ddof=1)
#     color = 'C{}'.format(i)
#     ax.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
#     ax.plot(dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
# ax.legend()
# ax.set_xlabel('number of dimensions')
# ax.set_ylabel('VAMP2 score')
# fig.tight_layout()
# plt.savefig(prefix+"sidechain_vamp2.png")
# plt.clf()
