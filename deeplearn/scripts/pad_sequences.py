import gzip
import numpy as np
import time
from scipy.stats import pearsonr, spearmanr, rankdata, ks_2samp
import matplotlib.pyplot as plt
import matplotlib

subsampleValues = [0.05, 0.1, 0.2, 0.4]
for prob in subsampleValues:
    fileName = 'sequencesJun26Classification' + str(prob) + 'NegativeSubsample.fa'
    print "on file " + fileName
    sequences = open('../features/' + fileName, 'r')
    paddedSequences = open('../features/padded' + fileName, 'w')

    print "padding..."

    PADDED_LENGTH = 2000
    i = 0
    for line in sequences:
        if i % 2 == 0:
            paddedSequences.write(line)
        else:
            seq = line.strip()
            length = len(seq)
            padLength = (PADDED_LENGTH - length) / 2
            newSeq = 'N' * padLength + seq + 'N' * padLength
            if len(newSeq) < PADDED_LENGTH:
                newSeq += 'N'
            paddedSequences.write(newSeq + '\n')
        i += 1
