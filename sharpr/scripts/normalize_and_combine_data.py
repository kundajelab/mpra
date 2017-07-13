import numpy as np
import os

def norm_signal(dna_in, rna_out):
    if isinstance(dna_in, list):
        dna_in = 

	return np.log2(rna_out / dna_in)

data_path = os.environ.get("SHARPR_DATA_DIR")
files = os.listdir(data_path)

rna_files = [rna_file for rna_file in files if rna_file.find("mRNA") != -1]
dna_files = [dna_file for dna_file in files if dna_file.find("Plasmid") != -1]

dna_reads = {}
dna_reads['1'] = {'minp': [], 'sv40p': []}
dna_reads['2'] = {'minp': [], 'sv40p': []}

# verbosely creating lists of 
rna_reads = {}
rna_reads['1'] = {}
rna_reads['2'] = {}
rna_reads['1']['k562'] = {}
rna_reads['1']['hepg2'] = {}
rna_reads['2']['k562'] = {}
rna_reads['2']['hepg2'] = {}
rna_reads['1']['k562']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['1']['k562']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['1']['hepg2']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['1']['hepg2']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['2']['k562']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['2']['k562']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['2']['hepg2']['minp'] = {'rep1': [], 'rep2': []}
rna_reads['2']['hepg2']['minp'] = {'rep1': [], 'rep2': []}


