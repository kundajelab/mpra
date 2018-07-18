import numpy as np
import os
import string
nucleotideMap = string.maketrans('ACGT', 'TGCA')

def rev(seq):
    return seq[::-1]

def comp(seq):
    return seq.translate(nucleotideMap)

def rev_comp(seq):
    return rev(comp(seq))

val_chrs = ['chr8']
test_chrs = ['chr18']
chrs = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
train_chrs = chrs
train_chrs = [chrom for chrom in train_chrs if chrom not in val_chrs and chrom not in test_chrs]

functional = True
run_name = 'counts_poisson_sep7'
deeplearn_dir = os.environ.get("DL")
split_path = deeplearn_dir + "/splits/" + run_name + "/"
os.system("mkdir " + split_path)

dataMatrix = open(os.environ.get("SHARPR") + "/data/processed_data/sharprFullDataMatrixZNormedLfcWCounts.tsv").readlines()[1:]
#dataMatrix.readline()

trainSplit = open(split_path + "train_split.txt", 'w')
valSplit = open(split_path + "val_split.txt", 'w')
testSplit = open(split_path + "test_split.txt", 'w')

if functional:
    funcLabelsPath = "%s/labels/functional_%s/" % (deeplearn_dir, run_name)
    os.system("mkdir %s" % funcLabelsPath)
    labels = open(funcLabelsPath + 'labelMatrix.txt', 'w')
else:
    labels = open(deeplearn_dir + "/labels/labels_" + run_name + ".txt", 'w')
labels.write("name\t" +
             #  "k562_minp_rep1_count\tk562_minp_rep2_count\tk562_minp_avg_count\t" +
             #  "k562_sv40p_rep1_count\tk562_sv40p_rep2_count\tk562_sv40p_avg_count\t" +
             #  "hepg2_minp_rep1_count\thepg2_minp_rep2_count\thepg2_minp_avg_count\t" +
             #  "hepg2_sv40p_rep1_count\thepg2_sv40p_rep2_count\thepg2_sv40p_avg_count\t" +
             #  "k562_minp_rep1_norm\tk562_minp_rep2_norm\tk562_minp_avg_norm\t" +
             #  "k562_sv40p_rep1_norm\tk562_sv40p_rep2_norm\tk562_sv40p_avg_norm\t" +
             #  "hepg2_minp_rep1_norm\thepg2_minp_rep2_norm\thepg2_minp_avg_norm\t" +
             #  "hepg2_sv40p_rep1_norm\thepg2_sv40p_rep2_norm\thepg2_sv40p_avg_norm\n"
             "k562_minp_rep1\tk562_minp_rep2\tk562_minp_avg\t" +
             "k562_sv40p_rep1\tk562_sv40p_rep2\tk562_sv40p_avg\t" +
             "hepg2_minp_rep1\thepg2_minp_rep2\thepg2_minp_avg\t" +
             "hepg2_sv40p_rep1\thepg2_sv40p_rep2\thepg2_sv40p_avg\n")

#  weights = open(deeplearn_dir + "/weights/weights_" + run_name + ".txt", 'w')
#  weights.write("name\t" +
              #  "k562_minp_rep1\tk562_minp_rep2\tk562_minp_avg\t" +
              #  "k562_sv40p_rep1\tk562_sv40p_rep2\tk562_sv40p_avg\t" +
              #  "hepg2_minp_rep1\thepg2_minp_rep2\thepg2_minp_avg\t" +
              #  "hepg2_sv40p_rep1\thepg2_sv40p_rep2\thepg2_sv40p_avg\n")
              
features = open(deeplearn_dir + "/features/sequences_" + run_name + ".fa", 'w')

counts = [0, 0, 0, 0]
num_fragments = len(dataMatrix)
print_levels = np.arange(0, num_fragments, num_fragments / 10, dtype = np.uint32)
for (i, line) in enumerate(dataMatrix):
    if i in print_levels:
        print("On fragment " + str(i) + " / " + str(num_fragments))
    line = line.strip().split('\t')
    
    k562_minp_rep1 = float(line[19])
    k562_minp_rep2 = float(line[20])
    # check if the datapoint is reproducible OR large fold-change
    useData = (abs(k562_minp_rep1 - k562_minp_rep2) < 0.5 or abs(k562_minp_rep1 + k562_minp_rep2)/2.0 > 1.7)

    for ori in ('n', 'rc'):
        chrom = line[1]
        if chrom in chrs:
            counts[3] += 1
    
        fragmentName = line[0] + '_' + ori
        if chrom in train_chrs:
            trainSplit.write(fragmentName + '\n')
            counts[0] += 1
        if chrom in val_chrs and ori == 'n' and (useData or not useData): # for val, don't get the revcomp sequence too
            valSplit.write(fragmentName + '\n')
            counts[1] += 1
        if chrom in test_chrs and ori in 'n':
            testSplit.write(fragmentName + '\n')
            counts[2] += 1
    
        if chrom in train_chrs or chrom in val_chrs or chrom in test_chrs:
            #  labels.write(fragmentName + '\t' + '\t'.join(line[19:31]) + '\n')
            labels.write(fragmentName + '\t' + '\t'.join(line[7:19]) + '\n')
            
            # assign weights to each sample's avg data, according to sigmoid(1 / (rep1 - rep2))
            #  from scipy.stats import logistic
            #  lbls = np.array(line[19:31]).astype(np.float)
            #  example_weights = np.ones(len(lbls))
            #  for idx in range(2, len(example_weights), 3):
                #  example_weights[idx] = logistic.cdf(np.reciprocal(np.abs(lbls[idx-2] - lbls[idx-1])))
            #  weights.write(fragmentName + '\t' + '\t'.join(example_weights.astype("string")) + '\n')

            seq = line[3]
            features.write('>' + fragmentName + '\n')
            if ori == 'rc':
                features.write(rev_comp(line[3]) + '\n')    
            else:
                features.write(line[3] + '\n')

trainSplit.close()
valSplit.close()
testSplit.close()
labels.close()
#  weights.close()
features.close()

print "train size = " + str(counts[0]) + " fragments"
print "val size = " + str(counts[1]) + " fragments"
print "test size = " + str(counts[2]) + " fragments"
print "total data size = " + str(counts[3]) + " fragments"
