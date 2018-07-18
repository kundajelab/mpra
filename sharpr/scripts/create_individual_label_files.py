import numpy as np
import pandas as pd
import os

base_path = os.environ.get("DL") + '/labels/functional_counts_poisson_sep7/'
data_path = base_path + 'labelMatrix.txt'
#  base_path = os.environ.get("DL") + '/weights/upweightends_aug10/'
#  data_path = base_path + 'replicatequality_only.txt'
allLabels = pd.read_csv(data_path,
                        dtype = 'string',
                        delimiter = '\t',
                        header = None,
                        #  max_rows = 10
                       ).values

for i in range(allLabels.shape[1] - 1):
    label_name = allLabels[0, i+1]
    np.savetxt(fname = base_path + label_name + '.txt',
               X = allLabels[1:, [0,i+1]],
               fmt = '%s',
               delimiter = '\t',
               header = 'name\t' + label_name,
               comments = ''
              )
