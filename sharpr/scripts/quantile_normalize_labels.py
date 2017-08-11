import os
import numpy as np
from scipy.stats import rankdata, spearmanr

def quantile_normalize(data):
    data_ranked = np.array([(rankdata(arr, method = 'min') - 1) for arr in data])
    data_sorted = np.array([np.sort(arr) for arr in data])
    average_by_rank = np.mean(data_sorted, axis = 0)
    data_normed = np.array([np.take(average_by_rank, ranks) for ranks in data_ranked])
    return data_normed

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
dataMatrixQNormed = dataMatrix
print dataMatrixQNormed.shape

expNames = ['k562_minp', 'hepg2_minp', 'k562_sv40p', 'hepg2_sv40p']

for (i, exp) in enumerate(expNames):
    exp_index = cols_to_indices[exp + '_rep1_count']
    rep1 = dataMatrix[:, exp_index].astype(np.float)
    rep2 = dataMatrix[:, exp_index + 1].astype(np.float)
    promoter = exp.split('_')[1]
    dna = dataMatrix[:, cols_to_indices['dna_' + promoter + '_count']].astype(np.float)
    quantNormedLfc = quantile_normalize([np.log2(rep1 + 1) - np.log2(dna + 1), 
                                         np.log2(rep2 + 1) - np.log2(dna + 1)])
    print(exp + " spearmanr: " +  str(spearmanr(quantNormedLfc[0], quantNormedLfc[1]))) # sanity check
    dataMatrixQNormed[:, exp_index + 12] = quantNormedLfc[0]
    dataMatrixQNormed[:, exp_index + 1 + 12] = quantNormedLfc[1]
    dataMatrixQNormed[:, exp_index + 2 + 12] = np.mean(quantNormedLfc, axis = 0)

np.savetxt(fname = base_path + "sharprFullDataMatrixQNormed.tsv",
           X = dataMatrixQNormed,
           fmt = '%s',
           delimiter = '\t',
           header = header,
           comments = '')
