
import sys
import pysam
from Bio import SeqIO
from pyfaidx import Fasta
from collections import Counter

comp_nts = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def rev_comp(seq):
	return ''.join([comp_nts[seq[i]] for i in range(len(seq)-1,-1,-1)])

def filter_fivep_reads(fivep_to_reads, read_file, out_path, info_out):

	read_sites = {}
	read_sequences = {}
	site_coords = {}
	read_fasta, read_sequences = Fasta(read_file, as_raw=True), {}
	in_bam = pysam.AlignmentFile(fivep_to_reads, 'rb')
	for read in in_bam.fetch():
		num_mismatch, read_cig = read.get_tag('XM'), read.cigarstring
		if num_mismatch < 1 and 'I' not in read_cig and 'D' not in read_cig:
			fivep_site, rid = read.query_name, read.reference_name
			fivep_site = fivep_site[:-3]
			if rid not in read_sites:
				read_sites[rid] = set()
				site_coords[rid] = {}
			read_sites[rid].add(fivep_site)
			read_sequences[rid] = read_fasta[rid][:]
			site_coords[rid][fivep_site] = (read.reference_start, read.reference_end, read.is_reverse)
	in_bam.close()

	with open(out_path, 'w') as out_file:
		with open(info_out, 'w') as out_info:
			out_info.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tfivep_first\tread_fivep_start\tread_fivep_end\n')
			for rid in read_sites:
				fivep_info = [site_coords[rid][fp] for fp in site_coords[rid]]
				if len(set(fivep_info)) == 1:
					fivep_start, fivep_end, is_reverse = fivep_info[0]
					read_seq = read_sequences[rid]
					if not is_reverse:
						trim_seq, fivep_seq = read_seq[:fivep_start], read_seq[fivep_start:fivep_end]
					else:
						trim_seq, fivep_seq = read_seq[fivep_end:], rev_comp(read_seq[fivep_start:fivep_end])
					if len(trim_seq) >= 20:
						out_file.write('>{}\n{}\n'.format(rid, trim_seq))
						fivep_sites = ','.join(site_coords[rid].keys())
						out_info.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(rid, read_seq, fivep_seq, fivep_sites, is_reverse, fivep_start, fivep_end))
					

if __name__ == '__main__' :

	read_file, fivep_to_reads, out_path, info_out = sys.argv[1:]
	filter_fivep_reads(fivep_to_reads, read_file, out_path, info_out)



			

