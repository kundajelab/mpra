import gzip
import numpy as np
import time
from scipy.stats import pearsonr, spearmanr, rankdata, ks_2samp
import matplotlib.pyplot as plt
import matplotlib
import os
import sys

negativeSubsamples = [0.05, 0.1, 0.2, 0.4]

# runName = 'Jun22ClassificationWithoutQc'
chrs = ['chr' + str(i) for i in range(1, 23)]
valChrs = ['chr8']
testChrs = ['chr18']
trainChrs = chrs
trainChrs = [chrom for chrom in trainChrs if chrom not in valChrs and chrom not in testChrs]
# print trainChrs

for i in range(len(negativeSubsamples)):
    prob = negativeSubsamples[i]
    runName = 'Jun26Classification' + str(prob) + 'NegativeSubsample'
    os.system("mkdir ../splits/" + runName + "/")
    dataMatrix = open('../data/classificationSignals' + str(prob) + 'NegativeSubsample.bed', 'r').readlines()[2:]
    # print len(dataMatrix)
    # sys.exit(0)
    # header = dataMatrix.readline()
    # spacer = dataMatrix.readline()

    trainSplit = open('../splits/' + runName + '/train_split.txt', 'w')
    valSplit = open('../splits/' + runName + '/val_split.txt', 'w')
    testSplit = open('../splits/' + runName + '/test_split.txt', 'w')

    # runName = 'Jun22ClassificationWithoutQc'
    classificationLabels = open('../labels/labels' + runName + '.txt', 'w')
    classificationLabels.write('name\tavg\trep1\trep2\n')
    # runName = 'Jun22RegressionWithoutQc'
#     regressionLabels = open('../labels/labels' + runName + '.txt', 'w')
#     regressionLabels.write('name\tavg\trep1\trep2\n')

    counts = [0,0,0,0]
    for line in dataMatrix:
        line = line.strip().split('\t')
        chrom = line[0]
        fragmentName = line[3]
        if chrom in chrs:
            counts[3] += 1

        # add fragment names to their respective split files
        if chrom in trainChrs:
            trainSplit.write(fragmentName + '\n')
            counts[0] += 1
        elif chrom in valChrs:
            valSplit.write(fragmentName + '\n')
            counts[1] += 1
        elif chrom in testChrs:
            testSplit.write(fragmentName + '\n')
            counts[2] += 1

        avgSignal = float(line[4])
        rep1Signal = float(line[6])
        rep2Signal = float(line[7])

        # if the fragment was in any of the relevant splits, then add it to the labels file
        if chrom in trainChrs or chrom in valChrs or chrom in testChrs:
            classificationLabel = fragmentName + '\t' + str(int(avgSignal > 0)) + '\t' + \
                                  str(int(rep1Signal > 0)) + '\t' + str(int(rep2Signal > 0)) + '\n'
            classificationLabels.write(classificationLabel)

#             regressionLabel = fragmentName + '\t' + str(avgSignal) + '\t' + \
#                               str(rep1Signal) + '\t' + str(rep2Signal) + '\n'
#             regressionLabels.write(regressionLabel)

    trainSplit.close()
    valSplit.close()
    testSplit.close()
    classificationLabels.close()
#     regressionLabels.close()

    print "train size = " + str(counts[0]) + " fragments"
    print "validation size = " + str(counts[1]) + " fragments"
    print "test size = " + str(counts[2]) + " fragments"
    print "total data size = " + str(counts[3]) + " fragments"

