########################################
##2019-7-1
##chenzhi
########################################
FEATURECOUNTS="/home/chenzhi/biosoft/subread-1.6.4-Linux-x86_64/bin"
DIR='/home/chenzhi/fareeha'
REP_INDEX=["G1T1R1","G1T1R2","G1T1R3",
	   "G1T2R1","G1T2R2","G1T2R3",
           "G2T1R1","G2T1R2","G2T1R3",
           "G2T2R1","G2T2R2","G2T2R3",
           "G1T1L1","G1T1L2","G1T1L3",
           "G1T2L1","G1T2L2","G1T2L3",
           "G2T1L1","G2T1L2","G2T1L3",
           "G2T2L1","G2T2L2","G2T2L3"]
INDEX_HT2="/home/chenzhi/fareeha/ref/watermelon"
REF_Genome="/home/chenzhi/fareeha/ref/watermelon_97103_v2.genome.fa"
GTF = "/home/chenzhi/fareeha/ref/watermelon_v2.gtf"


rule all:
	input:
		expand("{dir}/alignment/{rep}_sorted.bam",rep=REP_INDEX,dir=DIR),
		expand("{dir}/count/{rep}.count.tsv",rep=REP_INDEX,dir=DIR)
		#expand("{dir}/alignment/{rep}.sam",rep=REP_INDEX,dir=DIR)	


#rule Unpack1:
#	input:		
#		"{dir}/rawdata/{rep}_R1.gz"		
#	
#       output:
#		"{dir}/data/{rep}_R1"
#
#
#	shell:
#		"gunzip {input} {output}"
#
#rule Unpack2:
#	input:
#		"{dir}/rawdata/{rep}_R2.gz"
#	
#	output:
#		"{dir}/data/{rep}_R2"
#	
#	shell:
#		"gunzip {input} >{output}"
#
rule Alignment:
	input: 
		"{dir}/rawdata/{rep}_R1",
		"{dir}/rawdata/{rep}_R2"

	output:
		"{dir}/alignment/{rep}.sam"

	log:
		"{dir}/alignment/log/{rep}.alignment.log"
	
	threads:16

	shell:
		"hisat2 -p {threads} -x {INDEX_HT2} -1 {input[0]} -2 {input[1]} -S {output} 2>{log}"


rule sambam:	
	input:	
		"{dir}/alignment/{rep}.sam"

	output:
		"{dir}/alignment/{rep}_sorted.bam"
	
	threads:16

	shell:
		"samtools view -bS {input} |"
		" samtools sort -@ {threads} - -o {output}"



rule counts:
	input:
		"{dir}/alignment/{rep}_sorted.bam"

	output:
		"{dir}/count/{rep}.count.tsv"

	log:
		"{dir}/count/log/{rep}.summary.log"

	threads: 16

	shell:
		"{FEATURECOUNTS}/featureCounts -T {threads} -p -t exon -g gene_id -a {GTF} -o {output} {input} 2>{log}"
		

#rule cutadapt:
#	input:
#			"raw.fastq/ATAC_seq_{rep}_R1.fq.gz",
#			"raw.fastq/ATAC_seq_{rep}_R2.fq.gz"
#	output:
#			"fix.fastq/ATAC_seq_{rep}_R1.fq.gz",
#			"fix.fastq/ATAC_seq_{rep}_R2.fq.gz"
#	log:
#			"fix.fastq/ATAC_seq_{rep}_cutadapt.log"
#	shell:
#			"cutadapt  -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1"
#
#
#
#rule bt2_mapping: 
#	input:
#			"fix.fastq/ATAC_seq_{rep}_R1.fq.gz",
#			"fix.fastq/ATAC_seq_{rep}_R2.fq.gz"
#	output:
#			"bam/ATAC_seq_{rep}_bt2_hg38.sam"
#	log:
#			"bam/ATAC_seq_{rep}_bt2_hg38.log"
#	shell:
#			"bowtie2 -x {INDEX_BT2} -p 20 -1 {input[0]} -2 {input[1]} -s {output} > {log} 2>&1"
#
#
#rule bam_file_sort:
#	input:
#			"bam/ATAC_seq_{rep}_bt2_hg38.sam"
#	output:
#			"bam/ATAC_seq_{rep}_bt2_hg38_sort.bam"
#	log:
#			"bam/ATAC_seq_{rep}_bt2_hg38.sort.log"
#	shell:
#			"samtools sort -O BAM -o {output} -T {output}.temp -@ 4 -m 2G {input}"
#
#rule remove_duplication:
#	input:
#			"bam/ATAC_seq_{rep}_bt2_hg38_sort.bam"
#	output:
#			"bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.bam",
#			"bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.matrix"
#	log:
#			"bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.log"
#	shell:
#			"java -Xms5g -Xms5g -XX:ParallelGCThreads=4 \
#			- jar {PICARD} MarkDuplicates \
#			I={input} O={output[0]} M={output[1]} \
#			ASO=coordinate RMEOVE_DUPLICATES=true 2>{log}"
#
#rule call_peak:
#	input:
#			"bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.bam"
#	output:
#			"macs2_result/ATAC_seq_{rep}_peaks.narrowPeak"
#	params:
#			"ATAC_seq_{rep}",
#			"macs2_result"
#	log:
#			"macs2_result/ATAC_seq_{rep}_peaks.narrowPeak.log"
#	shell:
#			"macs2 callpeak -t {input} -f BAM -g hs --outdir {params[1]} -n {params[0]} -m 2 100 >{log} 2>&1"
#
