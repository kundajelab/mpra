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
