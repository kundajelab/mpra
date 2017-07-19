import numpy as np
import os

val_chrs = ['chr8']
test_chrs = ['chr18']
chrs = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
train_chrs = chrs
train_chrs = [chrom for chrom in train_chrs if chrom not in val_chrs and chrom not in test_chrs]

run_name = 'sharpr_regression_znormed_jul18'
deeplearn_dir = os.environ.get("DL")
split_path = deeplearn_dir + "/splits/" + run_name + "/"
os.system("mkdir " + split_path)

dataMatrix = open(os.environ.get("SHARPR") + "/data/processed_data/sharprFullDataMatrixZNormed.tsv").readlines()[1:]
#dataMatrix.readline()

trainSplit = open(split_path + "train_split.txt", 'w')
valSplit = open(split_path + "val_split.txt", 'w')
testSplit = open(split_path + "test_split.txt", 'w')

labels = open(deeplearn_dir + "/labels/labels_" + run_name + ".txt", 'w')
labels.write("name\t" +
             "k562_minp_rep1_count\tk562_minp_rep2_count\tk562_minp_avg_count\t" +
             "k562_sv40p_rep1_count\tk562_sv40p_rep2_count\tk562_sv40p_avg_count\t" +
             "hepg2_minp_rep1_count\thepg2_minp_rep2_count\thepg2_minp_avg_count\t" +
             "hepg2_sv40p_rep1_count\thepg2_sv40p_rep2_count\thepg2_sv40p_avg_count\t" +
             "k562_minp_rep1_norm\tk562_minp_rep2_norm\tk562_minp_avg_norm\t" +
             "k562_sv40p_rep1_norm\tk562_sv40p_rep2_norm\tk562_sv40p_avg_norm\t" +
             "hepg2_minp_rep1_norm\thepg2_minp_rep2_norm\thepg2_minp_avg_norm\t" +
             "hepg2_sv40p_rep1_norm\thepg2_sv40p_rep2_norm\thepg2_sv40p_avg_norm\n")
              
features = open(deeplearn_dir + "/features/sequences_" + run_name + ".fa", 'w')

counts = [0, 0, 0, 0]
num_fragments = len(dataMatrix)
print_levels = np.arange(0, num_fragments, num_fragments / 10, dtype = np.uint32)
for (i, line) in enumerate(dataMatrix):
    if i in print_levels:
        print("On fragment " + str(i) + " / " + str(num_fragments))
    line = line.strip().split('\t')
    chrom = line[1]
    if chrom in chrs:
        counts[3] += 1
    
    fragmentName = line[0]
    if chrom in train_chrs:
        trainSplit.write(fragmentName + '\n')
        counts[0] += 1
    if chrom in val_chrs:
        valSplit.write(fragmentName + '\n')
        counts[1] += 1
    if chrom in test_chrs:
        testSplit.write(fragmentName + '\n')
        counts[2] += 1
    
    if chrom in train_chrs or chrom in val_chrs or chrom in test_chrs:
        labels.write(fragmentName + '\t' + '\t'.join(line[7:31]) + '\n')
        seq = line[3]
        features.write('>' + fragmentName + '\n')
        features.write(line[3] + '\n')

trainSplit.close()
valSplit.close()
testSplit.close()
labels.close()
features.close()

print "train size = " + str(counts[0]) + " fragments"
print "val size = " + str(counts[1]) + " fragments"
print "test size = " + str(counts[2]) + " fragments"
print "total data size = " + str(counts[3]) + " fragments"
