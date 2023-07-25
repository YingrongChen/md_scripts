import matplotlib.pyplot as plt
import numpy as np
import pyemma
import sys

topfile=('backbone_nontcr.pdb')
train_files=('backbone_nontcr_prod.nc')
ahelixCA=[202,206,210,214,218,222,226,230,234,238,242,246,250,254,258,262,266,270,274,278,282,286,290,294,298,302,306]
bhelixCA=[554,558,562,566,570,574,578,582,586,590,594,598,602,606,610,614,618,622,626,630,634,638,642,646,650,654,658,662,666,670,674,678,682,686]
peptideCA=[718,722,726,730,734,738,742,746,750,754,758,762,766,770,774]
feat = pyemma.coordinates.featurizer(topfile)
feat.add_backbone_torsions(selstr="")
feat.add_distances(ahelixCA, peptideCA)
feat.add_distances(bhelixCA, peptideCA)
print('number of features ', feat.dimension())
data = pyemma.coordinates.load(train_files, features=feat)
tica = pyemma.coordinates.tica(data)
cluster = pyemma.coordinates.cluster_kmeans(tica, k=100, max_iter=50)
its = pyemma.msm.its(cluster.dtrajs, lags=[1, 2, 3, 5, 7, 10], nits=3, errors='bayes')
pyemma.plots.plot_implied_timescales(its, ylog=False);
plt.savefig('msm_its.png')