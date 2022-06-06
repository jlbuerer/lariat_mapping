
import sys
import gzip
import pysam
from intervaltree import Interval, IntervalTree
from collections import Counter
from subprocess import run
from random import sample
from os.path import join

chroms = set(['chr'+str(i) for i in range(1,23)] + ['chrY', 'chrX', 'chrM'])
comp_nts = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def rev_comp(seq):
	return ''.join([comp_nts[seq[i]] for i in range(len(seq)-1,-1,-1)])

def parse_cigar(cig_str):
	cig_tuples, curr_str = [], ''
	for c in cig_str:
		if not c.isalpha():
			curr_str = curr_str + c
		else:
			cig_tuples.append((int(curr_str), c))
			curr_str = ''
	return cig_tuples

def parse_read_info(read_info_table):

	read_info = {}
	with open(read_info_table) as in_file:
		for line in in_file:
			items = line.strip().split('\t')
			if items[0] != 'read_id':
				read_info[items[0]] = items[1:]

	return read_info

def parse_gene_ranges(anno_file):

	if anno_file[-2:] == 'gz':
		in_file = gzip.open(anno_file, 'rt')
	else:
		in_file = open(anno_file)

	gene_ranges = {}
	for line in in_file:
		if line[0] != '#':
			chrom, _, feat, start, end, _, strand, _, info = line.strip().split('\t')
			if feat == 'gene':
				info = {e.split(' ')[0]:e.split(' ')[1].replace('\"', '') for e in info[:-1].split('; ')}
				if chrom not in gene_ranges:
					gene_ranges[chrom] = {s:IntervalTree() for s in ('-', '+')}
				gene_ranges[chrom][strand].add(Interval(int(start)-1, int(end), {'gene_name':info['gene_name']}))
	in_file.close()

	return gene_ranges

def parse_threep_lens(threep_lens):

	threep_lengths = {}
	with open(threep_lens) as in_file:
		for line in in_file:
			chrom, threep_site, strand, length = line.strip().split('\t')
			threep_lengths[(chrom, strand, int(threep_site))] = int(length)

	return threep_lengths

def filter_threep_reads(reads_to_threep, threep_lengths, read_info, gene_ranges, genome_file, out_dir, sample_name):

	threep_info = {}
	in_bam = pysam.AlignmentFile(reads_to_threep, 'rb')
	fail_reason = Counter()
	for read in in_bam.fetch():
		rid = read.query_name
		num_mismatch, read_cig = read.get_tag('XM'), parse_cigar(read.cigarstring)
		read_len = float(sum(c[0] for c in read_cig if c[1] in ('S', 'I', 'M')))
		if num_mismatch <= 5 and num_mismatch/read_len <= 0.1:
			indels = [c[0] for c in read_cig if c[1] in ('D', 'I')]
			if len(indels) == 0 or (len(indels) == 1 and indels[0] <= 3):
				threep_site = read.reference_name[:-3]
				chrom, start, end, strand = threep_site.split(';')
				if strand == '+':
					start = end
				if rid not in threep_info:
					threep_info[rid] = []
				read_seq = read.get_forward_sequence() if not read.is_reverse else rev_comp(read.get_forward_sequence())
				threep_info[rid].append((chrom, strand, start, read.reference_start, read.reference_end, read.is_reverse, read_seq, read.query_alignment_end-1, read.mapping_quality))

	alignment_types = Counter()
	potential_alignments = {}
	for rid in threep_info:
		max_score = max(s[-1] for s in threep_info[rid])
		top_alignments = [s for s in threep_info[rid] if s[-1]==max_score]
		read_seq, fivep_seq, fivep_sites, read_is_reverse, fivep_start, fivep_end = read_info[rid]
		read_is_reverse = True if read_is_reverse=='True' else False

		fivep_dict = {c:{s:set() for s in ['+', '-']} for c in chroms}
		for fp in fivep_sites.split(','):
			chrom, start, end, strand = fp.split(';')
			if strand == '+':
				fivep_dict[chrom][strand].add(int(start))
			else:
				fivep_dict[chrom][strand].add(int(end))
		
		top_alignments_filtered = []
		for align_info in top_alignments:
			threep_chrom, threep_strand, threep_site, threep_start, threep_end, threep_is_reverse, threep_read, bp_read_site, _ = align_info
			threep_site, threep_end, bp_read_site = int(threep_site), int(threep_end), int(bp_read_site)
			if len(fivep_dict[threep_chrom][threep_strand]) > 0:
				threep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(threep_site, threep_site+1)]
				same_gene_fivep = []
				for fp in fivep_dict[threep_chrom][threep_strand]:
					fivep_genes = [g.data['gene_name'] for g in gene_ranges[threep_chrom][threep_strand].overlap(fp, fp+1)]
					gene_matches = sum(1 for g in fivep_genes if g in threep_genes) > 0
					if gene_matches:
						same_gene_fivep.append(fp)

				if len(same_gene_fivep) == 1:
					fivep_site = same_gene_fivep[0]
					if (strand == '+' and fivep_site < threep_site) or (strand == '-' and fivep_site > threep_site):
						if read_is_reverse == threep_is_reverse:
							if threep_strand == '+':
								bp_site = threep_end + threep_site-threep_lengths[(threep_chrom, threep_strand, threep_site)]-1
							else:
								bp_site = (threep_lengths[(threep_chrom, threep_strand, threep_site)]-threep_end) + threep_site
							if (strand == '+' and fivep_site < bp_site) or (strand == '-' and fivep_site > bp_site):
								bp_nt = threep_read[bp_read_site]
								if rid not in potential_alignments:
									potential_alignments[rid] = []
								potential_alignments[rid].append([read_seq, threep_chrom, threep_strand, fivep_site, read_is_reverse, fivep_start, fivep_end, threep_site, bp_site, bp_nt])

	temp_bp_bed, temp_bp_seq = join(out_dir, 'temp_bp_seqs.bed'), join(out_dir, 'temp_bp_seqs.txt')
	with open(join(out_dir, sample_name+'_final_info_table.txt'), 'w') as out_file:
		with open(join(out_dir, sample_name+'_lariat_data_table.txt'), 'w') as lar_out:
			out_file.write('read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tfivep_read_start\tfivep_read_end\t')
			out_file.write('threep_site\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window\n')
			for rid in potential_alignments:
				align_mismatch = {True:[], False:[]}
				for align_info in potential_alignments[rid]:
					read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_start, fivep_end, threep_site, bp_site, read_bp_nt = align_info
					temp_file = open(temp_bp_bed, 'w')
					if strand == '+':
						bp_start, bp_end = bp_site-4, bp_site+6
					else:
						bp_start, bp_end = bp_site-5, bp_site+5
					temp_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrom, bp_start, bp_end, '{};{};{}'.format(chrom, bp_site, strand), 0, strand))
					temp_file.close()
					run('bedtools getfasta -fi {} -bed {} -fo {} -nameOnly -s -tab'.format(genome_file, temp_bp_bed, temp_bp_seq).split(' '))
					temp_file = open(temp_bp_seq)
					name, genomic_bp_window = temp_file.readline().strip().split()
					temp_file.close()
					genomic_bp_window = genomic_bp_window.upper()
					genomic_bp_nt = genomic_bp_window[4]
					align_mismatch[genomic_bp_nt!=read_bp_nt].append(align_info+[genomic_bp_nt, genomic_bp_window])
					
				if len(align_mismatch[True]) > 0:
					output = [rid] + sample(align_mismatch[True], 1)[0]
				else:
					output = [rid] + sample(align_mismatch[False], 1)[0]
				out_file.write('\t'.join([str(e) for e in output]) + '\n')
				
				read_seq, chrom, strand, fivep_site, read_is_reverse, fivep_start, fivep_end, threep_site, bp_site, _, _, genomic_bp_window = output[1:]
				if read_is_reverse:
					read_seq = rev_comp(read_seq)
				output = [sample_name, 'exact', rid, read_seq, chrom, strand, fivep_site, threep_site, bp_site, genomic_bp_window]
				lar_out.write('\t'.join([str(e) for e in output]) + '\n')

	run('rm {} {}'.format(temp_bp_bed, temp_bp_seq).split(' '))


if __name__ == '__main__' :

	reads_to_threep, threep_lens, fivep_info_table, anno_file, genome_file, out_dir, sample_name = sys.argv[1:]
	gene_ranges, read_info = parse_gene_ranges(anno_file), parse_read_info(fivep_info_table)
	threep_lengths = parse_threep_lens(threep_lens)
	filter_threep_reads(reads_to_threep, threep_lengths, read_info, gene_ranges, genome_file, out_dir, sample_name)




