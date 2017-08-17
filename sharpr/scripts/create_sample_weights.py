import os
import numpy as np
from scipy.stats import logistic

base_path = os.environ.get("SHARPR") + "/data/processed_data/"
data_path = base_path + "sharprFullDataMatrixLfc.tsv"
dataMatrix = np.genfromtxt(fname = data_path,
                           dtype = 'string',
                           delimiter = '\t',
                           skip_header = 0,
                           #  max_rows = 10
                           )
header = dataMatrix[0][0] + '\t' + '\t'.join([task.replace("_norm", "") for task in dataMatrix[0][19:31]])
weightMatrix = np.ndarray(shape = (len(dataMatrix[1:]), 13), dtype = 'object')
weightMatrix[:, 0] = dataMatrix[1:, 0]
weightMatrix[:, 1:] = np.ones(shape=(len(weightMatrix), 12)) / 2
for idx in range(3, weightMatrix.shape[1], 3):
    # Downweight examples that are inconsistent across replicates
    weightMatrix[:, idx] = logistic.cdf(
                           np.reciprocal(
                           np.abs(
                           (dataMatrix[1:, 18+idx-2].astype(np.float) - dataMatrix[1:, 18+idx-1].astype(np.float) + 1e-6) # / (dataMatrix[1:, 18+idx-2].astype(np.float) + 1e-6)
                                  )))
    # Since sigmoid(x) is from 0.5 to 1 for x > 0, scale the range to 0 to 1
    weightMatrix[:, idx] = 1 - 2*(1 - weightMatrix[:, idx].astype(np.float))
    # Upweight repressive fragments
    repressiveFragments = dataMatrix[1:, 18+idx].astype(np.float) < -1 
    weightMatrix[:, idx][repressiveFragments] = 1*weightMatrix[:, idx][repressiveFragments]
    # Upweight the more "extreme" fragments to prevent the model from just predicting the mean
    extremeFragments = np.abs(dataMatrix[1:, 18+idx].astype(np.float)) > 30
    weightMatrix[:, idx][extremeFragments] = weightMatrix[:, idx][extremeFragments] * np.square(dataMatrix[1:, 18+idx][extremeFragments].astype(np.float))
    # Set all middle examples to weight 0 just to sanity check that weighting is working
    #  middleFragments = np.abs(dataMatrix[1:, 18+idx].astype(np.float)) >= 0
    #  weightMatrix[:, idx][middleFragments] = np.zeros(len(weightMatrix[:, idx][middleFragments]))
    # Convert float array to string
    # weightMatrix[:, idx] = weightMatrix[:, idx].astype("string")
    # Force the replicates to use the same weights as the average
    weightMatrix[:, idx-2] = weightMatrix[:, idx]
    weightMatrix[:, idx-1] = weightMatrix[:, idx]

weightMatrix = weightMatrix.astype("string")

weightFile = open(os.environ.get("DL") + "/weights/upweightends_aug10/replicatequality_only.txt", 'w')
weightFile.write(header + '\n')

for fragment in weightMatrix:
    fragment_n = list(fragment)
    fragment_rc = list(fragment)
    fragment_n[0] += '_n'
    weightFile.write('\t'.join(fragment_n) + '\n')
    fragment_rc[0] += '_rc'
    weightFile.write('\t'.join(fragment_rc) + '\n')
weightFile.close()

#  np.savetxt(fname = os.environ.get("DL") + "/weights/upweightweights_3x_weighted_aug7_norc.txt",
           #  X = weightMatrix,
           #  fmt = '%s',
           #  delimiter = '\t',
           #  header = header,
           #  comments = '')

#  f = open(os.environ.get("DL") + "/weights/weights_3x_weighted_aug7_norc.txt")
#  with_rc = open(os.environ.get("DL") + "/weights/weights_3x_weighted_aug7.txt", 'w')
#  with_rc.write(f.readline())
#  for line in f:
    #  line = line.strip().split('\t')
    #  newLine = list(line)
    #  newLine2 = list(line)
    #  newLine[0] += '_n'
    #  with_rc.write('\t'.join(newLine) + '\n')
    #  newLine2[0] += '_rc'
    #  with_rc.write('\t'.join(newLine2) + '\n')
#  f.close()
#  with_rc.close()

#  os.system("rm %s/weights/weights_3x_weighted_aug7_norc.txt" % os.environ.get("DL"))
