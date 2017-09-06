import numpy as np
import gzip
import sys

def get_coords_from_name(name):
	name = name.strip()
	chrom = name.split(':')[0]
	try:
		start, end = name[name.find(':') + 1 :].split('-')[:2]
	except:
		print name
		print name[name.find(':') :]
		sys.exit("Address the error")
	if name.find('(') >= 0:
		end = end[:end.find('(')]
		otherInfo = name[name.find('(') + 1 : name.find(')')]
		# print chrom, start, end, otherInfo, name
		# sys.exit("Check correctness")
		return [chrom, start, end, otherInfo, name]
	else:
		return [chrom, start, end, '.', name]
		
	

def write_seq_names_to_bed(seqfile = None, bedfile = None):
	if seqfile == None or bedfile == None:
		sys.exit("You must specify both an input file of sequence names and an output file name for the resulting BED.")
		return
	bed_output_file = open(bedfile, 'w')
	if seqfile[-3:] == '.gz':
		seq_names = gzip.open(seqfile)
	else:
		seq_names = open(seqfile)
	for line in seq_names:
		chrom, start, end, strand, name = get_coords_from_name(line)
		strand = strand.split(',')[0]
		bed_output_file.write(chrom + '\t' + start + '\t' + end + '\t' + name + '\t' + "0" + '\t' + strand + '\n')
	bed_output_file.close()

# From https://stackoverflow.com/a/21739593
def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty                                                                        
    pvalues = array(pvalues) 
    n = float(pvalues.shape[0])                                                                           
    new_pvalues = empty(int(n))
    if correction_type == "Bonferroni":                                                                   
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":                                                            
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]                                                                                                                  
    return new_pvalues
