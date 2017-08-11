import os
import numpy as np
import scipy
from scipy.stats import rankdata

def gaussian_projection(x):
    percentiles = rankdata(x, "average") / len(x)
    max_idx = np.argmax(percentiles)
    print percentiles[max_idx]
    if percentiles[max_idx] == 1.0:
        percentiles[max_idx] = 1.0 - 1.0 / len(x+1)
    return scipy.stats.norm.ppf(percentiles)

base_path = os.environ.get("SHARPR") + "/data/processed_data/"
data_path = base_path + "sharprFullDataMatrixLfc.tsv"
dataMatrix = np.genfromtxt(fname = data_path,
                           dtype = 'string',
                           delimiter = '\t',
                           skip_header = 0)
header = '\t'.join(dataMatrix[0])
dataMatrix = dataMatrix[1:]

dataMatrixProjected = dataMatrix
for col in range(19, 31):
    # print "Performing Gaussian projection for %s" % header.split('\t')[col]
    dataMatrixProjected[:, col] = gaussian_projection(dataMatrixProjected[:, col].astype(np.float)).astype("string")

print dataMatrixProjected.shape
print np.mean(dataMatrixProjected[:, 7:31].astype(np.float))
print np.std(dataMatrixProjected[:, 7:31].astype(np.float))
print np.mean(dataMatrixProjected[:, 7].astype(np.float))
print np.std(dataMatrixProjected[:, 7].astype(np.float))
print np.mean(dataMatrixProjected[:, 30].astype(np.float))
print np.std(dataMatrixProjected[:, 30].astype(np.float))
np.savetxt(fname = base_path + "sharprFullDataMatrixLfcGaussianProjected.tsv",
           X = dataMatrixProjected,
           fmt = '%s',
           delimiter = '\t',
           header = header,
           comments = '')
