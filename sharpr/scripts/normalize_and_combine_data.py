import numpy as np
import os
from collections import OrderedDict 

def norm_signal(dna_in, rna_out):
    if isinstance(dna_in, list):
        dna_in = np.array(dna_in)
    if isinstance(rna_out, list):
        rna_out = np.array(rna_out)
	return np.log2(rna_out / dna_in)

data_path = os.environ.get("SHARPR_DATA_DIR")
raw_data_path = data_path + "/rawCounts/"
files = os.listdir(raw_data_path)

rna_files = [raw_data_path + rna_file for rna_file in files if rna_file.find("mRNA") != -1]
dna_files = [raw_data_path + dna_file for dna_file in files if dna_file.find("Plasmid") != -1]

dna_reads = {}
dna_reads['1'] = {'minp': {}, 'sv40p': {}}
dna_reads['2'] = {'minp': {}, 'sv40p': {}}

# verbosely creating lists of the DNA/RNA reads
rna_reads = {}
rna_reads['1'] = {}
rna_reads['2'] = {}
rna_reads['1']['k562'] = {}
rna_reads['1']['hepg2'] = {}
rna_reads['2']['k562'] = {}
rna_reads['2']['hepg2'] = {}
rna_reads['1']['k562']['minp'] = {'rep1': {}, 'rep2': {}}
rna_reads['1']['k562']['sv40p'] = {'rep1': {}, 'rep2': {}}
rna_reads['1']['hepg2']['minp'] = {'rep1': {}, 'rep2': {}}
rna_reads['1']['hepg2']['sv40p'] = {'rep1': {}, 'rep2': {}}
rna_reads['2']['k562']['minp'] = {'rep1': {}, 'rep2': {}}
rna_reads['2']['k562']['sv40p'] = {'rep1': {}, 'rep2': {}}
rna_reads['2']['hepg2']['minp'] = {'rep1': {}, 'rep2': {}}
rna_reads['2']['hepg2']['sv40p'] = {'rep1': {}, 'rep2': {}}

names_to_info = OrderedDict()

total_dna_rna_read_counts = {'dna': {'1': {'minp': 0, 'sv40p': 0},
                                     '2': {'minp': 0, 'sv40p': 0}},
                             'rna': {'1': {'k562': {'minp': {'rep1': 0, 'rep2': 0},
                                                    'sv40p': {'rep1': 0, 'rep2': 0},
                                                   },
                                           'hepg2': {'minp': {'rep1': 0, 'rep2': 0},
                                                     'sv40p': {'rep1': 0, 'rep2': 0},
                                                    },
                                          },
                                     '2': {'k562': {'minp': {'rep1': 0, 'rep2': 0},
                                                    'sv40p': {'rep1': 0, 'rep2': 0},
                                                    },                                                                                                                                                                           'hepg2': {'minp': {'rep1': 0, 'rep2': 0},                                                                                                                                                       'sv40p': {'rep1': 0, 'rep2': 0},                                                                                                                                                     }
                                          }
                                    }
                            }
                                       
PSEUDOCOUNT = 1.0

# Read in DNA counts w/o the +1 pseudocount, compute running sum of total number of reads
for dna_file in dna_files:
    print("Loading from file " + dna_file)
    f = open(dna_file)
    design = ('1' if dna_file.find("ScaleUpDesign1") >= 0 else '2')
    promoter = ("minp" if dna_file.find("minP") >= 0 else "sv40p")
    f.readline() # skip header
    for line in f:
        line = line.strip().split('\t')
        name = line[0]
        chrom_state = name.split('_')[1]
        chrom = name.split('_')[4]
        center_pos = name.split('_')[5]
        seq = line[1]
        barcode = line[2] # not used for anything, but storing it anyway
        read_count = float(line[3])
        if promoter == 'minp': # only need to create name to sequence dictionary for one set of experiments
            names_to_info[name] = [chrom, center_pos, seq, barcode, chrom_state, design]
        dna_reads[design][promoter][name] = read_count # + PSEUDOCOUNT
        total_dna_rna_read_counts['dna'][design][promoter] += read_count + PSEUDOCOUNT
    f.close()

# Read in RNA counts w/o the +1 pseudocount, compute running sums of total numbers of reads
for rna_file in rna_files:
    print("Loading from file " + rna_file)
    f = open(rna_file)
    design = ('1' if rna_file.find("ScaleUpDesign1") >= 0 else '2')
    celltype = ('k562' if rna_file.find("K562") >= 0 else 'hepg2')
    promoter = ('minp' if rna_file.find('minP') >= 0 else 'sv40p')
    replicate = ('rep1' if rna_file.find('Rep1') >= 0 else 'rep2')
    f.readline() # skip header
    for line in f:
        line = line.strip().split('\t')
        name = line[0]
        seq = line[1]
        barcode = line[2] # not used for anything, but storing it anyway
        read_count = float(line[3])
        rna_reads[design][celltype][promoter][replicate][name] = read_count # + PSEUDOCOUNT
        total_dna_rna_read_counts['rna'][design][celltype][promoter][replicate] += read_count + PSEUDOCOUNT
    f.close()

dataMatrix = open(data_path + "/processed_data/sharprFullDataMatrixLfcPooled.tsv", 'w')
print("Creating data matrix at path " + data_path + "/processed_data/sharprFullDataMatrixLfc.tsv")
dataMatrix.write("name\tchrom\tcenter_coord\tsequence\tbarcode\tchromatin_state\tdesign\t" +
                 "k562_minp_rep1_count\tk562_minp_rep2_count\tk562_minp_avg_count\t" + 
                 "k562_sv40p_rep1_count\tk562_sv40p_rep2_count\tk562_sv40p_avg_count\t" +
                 "hepg2_minp_rep1_count\thepg2_minp_rep2_count\thepg2_minp_avg_count\t" + 
                 "hepg2_sv40p_rep1_count\thepg2_sv40p_rep2_count\thepg2_sv40p_avg_count\t" +
                 "k562_minp_rep1_norm\tk562_minp_rep2_norm\tk562_minp_avg_norm\t" + 
                 "k562_sv40p_rep1_norm\tk562_sv40p_rep2_norm\tk562_sv40p_avg_norm\t" +
                 "hepg2_minp_rep1_norm\thepg2_minp_rep2_norm\thepg2_minp_avg_norm\t" + 
                 "hepg2_sv40p_rep1_norm\thepg2_sv40p_rep2_norm\thepg2_sv40p_avg_norm\t" + 
                 "dna_minp_count\tdna_sv40p_count\n")

# all signals are log2 normalized
num_fragments = len(names_to_info.keys())
print_levels = np.arange(0, num_fragments, num_fragments / 10, dtype = np.uint32)

import itertools
cell_types = ['k562','hepg2']
promoters = ['minp', 'sv40p']
reps = ['rep1', 'rep2', 'avg']
# PSEUDOCOUNT = 1

from scipy.stats import rankdata

#  def quantile_normalize(data):
    #  data_ranked = np.array([(rankdata(arr, method = 'min') - 1) for arr in data])
    #  data_sorted = np.array([np.sort(arr) for arr in data])
    #  average_by_rank = np.mean(data_sorted, axis = 0)
    #  data_normed = np.array([np.take(average_by_rank, ranks) for ranks in data_ranked])
    #  return data_normed

pooled = True

for (i, seq_name) in enumerate(names_to_info.keys()):
    if i in print_levels:
        print("On fragment " + str(i) + " / " + str(num_fragments))
    dataMatrix.write(seq_name+'\t')
    seq_info = names_to_info[seq_name]
    dataMatrix.write("\t".join(names_to_info[seq_name]) + '\t')
    design = seq_info[5]
    # rna_out and dna_in arrays follow the order of the signals in the header line
    seq_labels = np.zeros(24)
    j = 0
    for ct, p, rep in itertools.product(cell_types, promoters, reps):
        if rep != 'avg':
            rna_count = rna_reads[design][ct][p][rep][seq_name]
            dna_count = dna_reads[design][p][seq_name]
            seq_labels[j] = rna_count
            seq_labels[j+12] = np.log2(rna_count + 1) - np.log2(dna_count + 1)
            #  seq_labels[j+12] = (np.log2(rna_count + 1) - np.log2(dna_count + 1) - 
                                #  np.log2(total_dna_rna_read_counts['rna'][design][ct][p][rep]) + 
                                #  np.log2(total_dna_rna_read_counts['dna'][design][p]))
        if rep == 'avg':
            seq_labels[j] = (seq_labels[j-2] + seq_labels[j-1]) / 2.0
            if pooled:
                rna_count_rep1 = rna_reads[design][ct][p]['rep1'][seq_name]
                rna_count_rep2 = rna_reads[design][ct][p]['rep2'][seq_name]
                dna_count = dna_reads[design][p][seq_name]
                seq_labels[j+12] = (np.log2(rna_count_rep1 + rna_count_rep2 + 1) -
                                    np.log2(dna_count + 1))
            else:
                seq_labels[j+12] = (seq_labels[j+10] + seq_labels[j+11]) / 2.0
        j += 1
    dataMatrix.write('\t'.join(seq_labels.astype('string')) + '\t')
    dataMatrix.write(str(dna_reads[design]['minp'][seq_name]) + '\t' +
                     str(dna_reads[design]['sv40p'][seq_name]) + '\n')
#    rna_out = [(rna_reads[design]['k562']['minp']['rep1'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['k562']['minp']['rep2'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['k562']['sv40p']['rep1'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['k562']['sv40p']['rep2'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['hepg2']['minp']['rep1'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['hepg2']['minp']['rep2'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['hepg2']['sv40p']['rep1'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1'],
#               (rna_reads[design]['hepg2']['sv40p']['rep2'][seq_name]) / 
#               total_dna_rna_read_counts['rna'][design]['k562']['minp']['rep1']
#               ]
#    dna_in =  [(dna_reads[design]['minp'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['minp'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['sv40p'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['sv40p'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['minp'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['minp'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['sv40p'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp'],
#               (dna_reads[design]['sv40p'][seq_name]) / 
#               total_dna_rna_read_counts['dna'][design]['minp']
#               ]
#    norm_signals = norm_signal(dna_in, rna_out)           
#    dataMatrix.write('\t'.join(norm_signals.astype('string')) + '\n')
    #if i < 3:
    #    #print(seq_name+'\t'+"\t".join(names_to_info[seq_name]) + '\t'+'\t'.join(norm_signals.astype('string')) + '\n')
    #    print(dna_in)
    #    print(rna_out)
    #    print('\t'.join(norm_signals.astype('string')))
    #else:
    #    break
print("Finished creating data matrix")
dataMatrix.close()

        

