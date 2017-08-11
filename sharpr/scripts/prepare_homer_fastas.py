import numpy as np
import os

base_path = os.environ.get("SHARPR") + "/data/processed_data/"
data_path = base_path + "sharprFullDataMatrixLfcZNormed.tsv"
dataMatrix = np.genfromtxt(fname = data_path,
                           dtype = 'string',
                           delimiter = '\t',
                           skip_header = 0)
header = '\t'.join(dataMatrix[0])
header_names = header.strip().split('\t')
cols_to_indices = {col: header_names.index(col) for col in header_names}
dataMatrix = dataMatrix[1:]

name_col = cols_to_indices['name']
seq_col = cols_to_indices['sequence']
numPositives = 20000
numNegatives = 20000
numBackground = 50000
exps = ["K562_minP", "K562_SV40P", "HepG2_minP", "HepG2_SV40P"]
for exp in exps:
    print exp

    homerPositives = open(os.environ.get("SHARPR") + "/data/homer_files/" + exp + "_top%d_positives.fa" % numPositives, 'w')
    homerNegatives = open(os.environ.get("SHARPR") + "/data/homer_files/" + exp + "_bottom%d_negatives.fa" % numNegatives, 'w')
    homerBackground = open(os.environ.get("SHARPR") + "/data/homer_files/" + exp + "_%d_background.fa" % numBackground, 'w')

    exp_col_idx = cols_to_indices[exp.lower() + "_avg_norm"]
    signed_sorted = np.argsort(dataMatrix[:, exp_col_idx].astype(np.float))
    
    pos_idxs = signed_sorted[-1 * numPositives:]
    for idx in pos_idxs:
        homerPositives.write('>' + dataMatrix[idx, name_col] + '\n')
        homerPositives.write(dataMatrix[idx, seq_col] + '\n')
    homerPositives.close()    

    neg_idxs = signed_sorted[:numNegatives]
    for idx in neg_idxs:
        homerNegatives.write('>' + dataMatrix[idx, name_col] + '\n')
        homerNegatives.write(dataMatrix[idx, seq_col] + '\n')
    homerNegatives.close()

    unsigned_sorted = np.argsort(np.abs(dataMatrix[:, exp_col_idx].astype(np.float)))
    background_idxs = unsigned_sorted[:numBackground]
    for idx in background_idxs:
        homerBackground.write('>' + dataMatrix[idx, name_col] + '\n')
        homerBackground.write(dataMatrix[idx, seq_col] + '\n')
    homerBackground.close()
