import matplotlib.pyplot as plt
import numpy as np
import pyemma
import sys

topfile=('backbone_nontcr.pdb')
trajfile=('backbone_nontcr_prod_small.nc')
feat = pyemma.coordinates.featurizer(topfile)
ser=feat.select("(resSeq 189) and name CA")
notser=feat.select("(not resSeq 189) and name CA")
helixCA=feat.select("(resSeq 51 to 77 or resSeq 139 to 172) and name CA")
peptideCA=feat.select("resSeq 181 to 194 and name CA")
print(peptideCA)
feat.add_backbone_torsions()
print('number of backbone torsions features ', feat.dimension())
feat.add_distances(helixCA, indices2=peptideCA)
# feat.add_distances(helixCA)
print('number of features ', feat.describe())
print('number of features ', feat.dimension())
data = pyemma.coordinates.load(trajfile, features=feat)
tica = pyemma.coordinates.tica(data)
cluster = pyemma.coordinates.cluster_kmeans(tica, k=100, max_iter=50)
its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10, 13, 15, 17, 20], nits=3, errors='bayes')
pyemma.plots.plot_implied_timescales(its, outfile="msm_its.png", show_mle=False);
