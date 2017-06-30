import gzip
import numpy as np
import time
from scipy.stats import pearsonr, spearmanr, rankdata, ks_2samp
import matplotlib.pyplot as plt
import matplotlib

signals = open('../../sureseq/SUREfragmentsignal3.txt')
header = signals.readline()

negativeSubsamples = [0.05, 0.1, 0.2, 0.4]

rep1Data = []
rep2Data = []
rep1GoodData = [[] for p in range(len(negativeSubsamples))]
rep2GoodData = [[] for p in range(len(negativeSubsamples))]
goodSignals = [[] for p in range(len(negativeSubsamples))]

i = 0
for line in signals:
    if (i-1) % 1e7 == 0:
        print "On fragment %s / " %(i-1) + '153000000'
    line = line.strip().split('\t')
    rep1Signal = float(line[4])
    rep2Signal = float(line[5])
    rep1Data.append(rep1Signal)
    rep2Data.append(rep2Signal)
    seqLen = int(line[2]) - int(line[1])
    for j in range(len(negativeSubsamples)):
        prob = negativeSubsamples[j]
        if (seqLen <= 2000 and
           ((rep1Signal > 0 and rep2Signal > 0) or 
           (rep1Signal == 0 and rep2Signal == 0 and np.random.uniform() < prob))):
    #     if rep1Signal > 0 and rep2Signal > 0 and seqLen <= 2000:
            rep1GoodData[j].append(rep1Signal)
            rep2GoodData[j].append(rep2Signal)
            name = line[0] + ":" + line[1] + "-" + line[2] + '(' + line[3] + ')'
            avgSignal = float(line[6])
            goodSignals[j].append([line[0], line[1], line[2], name, avgSignal, line[3], rep1Signal, rep2Signal])
    i += 1

goodSignals = np.array(goodSignals)
print goodSignals.shape
numPositive = 2000000
for k in range(len(goodSignals)):
    prob = negativeSubsamples[k]
    print("Length of array with " + str(prob) + " subsampling = " + str(len(goodSignals[k])) + ". "
          "Expected length was approximately " + str(int(numPositive + prob * (i - numPositive))) + ".")

print goodSignals[0][:5]

header = 'chr\tstart\tend\tname\tavgSignal\tstrand\trep1Signal\trep2Signal\n'

for j in range(len(negativeSubsamples)):
    prob = negativeSubsamples[j]
    np.savetxt(fname = '../data/classificationSignals' + str(prob) + 'NegativeSubsample.bed',
               fmt = '%s',
               X = goodSignals[j],
               delimiter = '\t', 
               header = header,
               comments='# ')
