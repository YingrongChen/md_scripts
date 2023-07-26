import numpy as np
import pyemma
topfile=('backbone_nontcr.pdb')
feat = pyemma.coordinates.featurizer(topfile)
ahelix=feat.select("(resSeq 51 to 77) and name CA")
bhelix=feat.select("(resSeq 139 to 172) and name CA")
peptide=feat.select("resSeq 181 to 194 and name CA")
# bhelixCA=[546,550,554,558,562,566,570,574,578,582,586,590,594,598,602,606,610,614,618,622,626,630,634,638,642,646,650,654,658,662,666,670,674,678,682,686]
# ahelixCA=[202,206,210,214,218,222,226,230,234,238,242,246,250,254,258,262,266,270,274,278,282,286,290,294,298,302,306]
# peptideCA=[718,722,726,730,734,738,742,746,750,754,758,762,766,770,774]
# n = len(ahelixCA)
# matrix = np.zeros((3*n,2), dtype=int)
# for i in range(0,n):
#     matrix[i,:] = [bhelixCA[i], ahelixCA[n-1-i]]
# for i in range(0,n):
#     matrix[i+n,:] = [bhelixCA[i+1], ahelixCA[n-1-i]]
# for i in range(0,n):
#     matrix[i+2*n,:] = [bhelixCA[i+2], ahelixCA[n-1-i]]
# print(matrix)

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
print(matrix)