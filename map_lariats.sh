#!/bin/bash

#single read file (paired-end read files should be processed separately)
read_file=$1 
#all reads should be trimmed to a uniform length
read_len=$2
sample_name=$3
#bowtie genome index base name
bowtie_genome_index=$4
#Fasta files containing the 5' (first 20nts) and 3' (last 250nts) splice site regions of all human introns
fivep_fasta=$5
threep_fasta=$6
#File listing the length of the 3'SS regions for each intron
threep_lens=$7
#GTF annotation file containing the exon information for mapping genome
gtf_file=$8
#number of threads available
threads=$9
#output directory (this script outputs lariat mapping data to out_dir/lariat_reads)
out_dir="${10}"

if [[ ! -d $out_dir ]]; then
	mkdir $out_dir
fi

### Filter reads with >5% ambiguous characters
echo "$(date +'%m/%d %H:%M:%S') - Filtering reads with >5% ambiguous characters"
filtered_reads=$out_dir/filtered_reads.fq
filtered_reads=$read_file
max_Ns=$(echo "$read_len * 0.05" | bc)
max_Ns=$(printf "%.0f" $max_Ns)
bbduk.sh in=$read_file out=$filtered_reads maxns=$max_Ns t=$threads > /dev/null

### Map filtered reads to genome and keep unmapped reads
echo "$(date +'%m/%d %H:%M:%S') - Mapping reads and extracting unmapped reads"
unmapped_bam=$out_dir/unmapped_reads.bam
bowtie2 --threads $threads --end-to-end --sensitive --score-min L,0,-0.24 -k 1 --n-ceil L,0,0.05 -x $bowtie_genome_index -U $filtered_reads \
	| samtools view -b -h -f 4 - > $unmapped_bam

### Create fasta file of unmapped reads and build bowtie index of it
echo "$(date +'%m/%d %H:%M:%S') - Creating fasta file of unmapped reads" 
unmapped_fasta=$out_dir/unmapped_reads.fa
samtools fasta $unmapped_bam > $unmapped_fasta
samtools faidx $unmapped_fasta
echo "$(date +'%m/%d %H:%M:%S') - Building bowtie index of unmapped fasta"
bowtie2-build --threads $threads --large-index $unmapped_fasta $unmapped_fasta > /dev/null

### Make fasta file of all 5' splice sites (first intronic 20nts) and align them to unmapped reads
echo "$(date +'%m/%d %H:%M:%S') - Mapping 5' splice sites to reads"
fivep_to_reads=$out_dir/fivep_to_reads.bam
/datasets2/jbuerer/bin/bowtie2 --threads $threads --end-to-end --sensitive --no-unal -f --k 10000 -x $unmapped_fasta -U $fivep_fasta \
	| samtools view -b -h - > $fivep_to_reads

samtools sort -o $out_dir/temp.bam -@ $threads $fivep_to_reads
mv $out_dir/temp.bam $fivep_to_reads
samtools index -@ $threads $fivep_to_reads

### Extract reads with a mapped 5' splice site and trim it off
echo "$(date +'%m/%d %H:%M:%S') - Finding 5' read alignments and trimming reads"
fivep_info_table=$out_dir/fivep_info_table.txt
python filter_fivep_alignments.py $unmapped_fasta $fivep_to_reads $out_dir/fivep_mapped_reads_trimmed.fa $fivep_info_table

### Make fasta file of all 3' splice sites (last intronic 250nts) and build bowtie index of it
echo "$(date +'%m/%d %H:%M:%S') - Building bowtie index of 3' splice site fasta"
bowtie2-build --threads $threads $threep_fasta $threep_fasta > /dev/null

### Map 5p trimmed reads to 3p fasta
echo "$(date +'%m/%d %H:%M:%S') - Mapping 5' trimmed reads to 3p fasta file"
reads_to_threep=$out_dir/fivep_reads_trimmed_mapped_to_threep.bam
/datasets2/jbuerer/bin/bowtie2 --threads $threads --end-to-end --sensitive -k 10 --no-unal -f -x $threep_fasta -U $out_dir/fivep_mapped_reads_trimmed.fa \
	| samtools view -b -h - > $reads_to_threep

samtools sort -o $out_dir/temp.bam -@ $threads $reads_to_threep
mv $out_dir/temp.bam $reads_to_threep
samtools index -@ $threads $reads_to_threep

### Filter 3' splice site alignments and output info table with BP site
echo "$(date +'%m/%d %H:%M:%S') - Filtering 3' alignments and outputting final table"
python filter_threep_alignments.py $reads_to_threep $threep_lens $fivep_info_table $gtf_file $bowtie_genome_index".fa" $out_dir $sample_name

rm $filtered_reads $unmapped_bam $out_dir/unmapped_reads.fa* $out_dir/fivep_mapped_reads_trimmed.fa $reads_to_threep* $fivep_to_reads*


