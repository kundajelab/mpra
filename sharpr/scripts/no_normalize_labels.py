import os
import numpy as np
from scipy.stats import rankdata, spearmanr

base_path = os.environ.get("SHARPR") + "/data/processed_data/"
data_path = base_path + "sharprFullDataMatrix.tsv"
dataMatrix = np.genfromtxt(fname = data_path,
                           dtype = 'string',
                           delimiter = '\t',
                           skip_header = 0)
header = '\t'.join(dataMatrix[0])
header_names = header.strip().split('\t')
cols_to_indices = {col: header_names.index(col) for col in header_names}
dataMatrix = dataMatrix[1:]
dataMatrixLfc = dataMatrix
print dataMatrixLfc.shape

expNames = ['k562_minp', 'hepg2_minp', 'k562_sv40p', 'hepg2_sv40p']

for (i, exp) in enumerate(expNames):
    exp_index = cols_to_indices[exp + '_rep1_count']
    rep1 = dataMatrix[:, exp_index].astype(np.float)
    rep2 = dataMatrix[:, exp_index + 1].astype(np.float)
    promoter = exp.split('_')[1]
    dna = dataMatrix[:, cols_to_indices['dna_' + promoter + '_count']].astype(np.float)
    lfc = [np.log2(rep1 + 1) - np.log2(dna + 1), 
           np.log2(rep2 + 1) - np.log2(dna + 1),
           np.log2(rep1 + rep2 + 1) - np.log2(dna + 1)] # average is just pooled number of reads normed by input
    print(exp + " spearmanr: " +  str(spearmanr(lfc[0], lfc[1]))) # sanity check
    dataMatrixLfc[:, exp_index + 12] = lfc[0]
    dataMatrixLfc[:, exp_index + 1 + 12] = lfc[1]
    dataMatrixLfc[:, exp_index + 2 + 12] = lfc[2]

np.savetxt(fname = base_path + "sharprFullDataMatrixLfcAvgPooled.tsv",
           X = dataMatrixLfc,
           fmt = '%s',
           delimiter = '\t',
           header = header,
           comments = '')
