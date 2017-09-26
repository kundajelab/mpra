import os
import numpy as np
from scipy.stats import zscore

base_path = os.environ.get("SHARPR") + "/data/processed_data/"
data_path = base_path + "sharprFullDataMatrixLfcPooled.tsv"
dataMatrix = np.genfromtxt(fname = data_path,
                           dtype = 'string',
                           delimiter = '\t',
                           skip_header = 0)
header = '\t'.join(dataMatrix[0])
dataMatrix = dataMatrix[1:]
dataMatrixZNormed = dataMatrix
print dataMatrixZNormed.shape
dataMatrixZNormed[:, 19:] = np.transpose(zscore(np.transpose(dataMatrixZNormed[:, 19:].astype(np.float)), axis = 1))
print dataMatrixZNormed.shape
print np.mean(dataMatrixZNormed[:, 7:31].astype(np.float))
print np.std(dataMatrixZNormed[:, 7:31].astype(np.float))
print np.mean(dataMatrixZNormed[:, 7].astype(np.float))
print np.std(dataMatrixZNormed[:, 7].astype(np.float))
print np.mean(dataMatrixZNormed[:, 30].astype(np.float))
print np.std(dataMatrixZNormed[:, 30].astype(np.float))
np.savetxt(fname = base_path + "sharprFullDataMatrixZNormedLfcWCounts.tsv",
           X = dataMatrixZNormed,
           fmt = '%s',
           delimiter = '\t',
           header = header,
           comments = '')
