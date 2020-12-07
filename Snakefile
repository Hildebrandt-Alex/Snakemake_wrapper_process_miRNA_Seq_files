import os
import pandas as pd

configfile: 'config.yaml'

samples = pd.read_table("samples.csv", header=0, sep=',', index_col=0)
lengths = list(range(config['min_length'],config['max_length']+1))


hairpin = os.path.join(config['reference_folder'],config['hairpin_ref_file'])
mature = os.path.join(config['reference_folder'],config['mature_db_file'])


rule all:
	input: 
#  	NGS-Seq_Data_and_Output_I_output
		expand('Seq_Data_and_Output/{sample}_Trimmed.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_Untrimmed.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_Trimming_Output.txt', sample=samples.index),
		expand('Seq_Data_and_Output/FastQC/{sample}_Trimmed_fastqc.html', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_rRNA_unmapped.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_rRNA_readcount.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_r_tRNA_unmapped.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_tRNA_readcount.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_isomir_readcount.txt', sample=samples.index),
		'Seq_Data_and_Output/isomir_readcount_all.txt',
		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_snRNA_readcount.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_snoRNA_readcount.txt', sample=samples.index),
# 	plotting Information_output
		expand('Seq_Data_and_Output/{sample}_short_noadapter.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_trimmedreads_basic_analyses.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_trimmedreads_lengths.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_Short.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_tmp.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_short_lengths.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_rRNA_unmapped_basic_analyses.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_rRNA_unmapped_lengths.txt', sample=samples.index),
		expand('Seq_Data_and_Output/{sample}_r_tRNA_unmapped_basic_analyses.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_unmapped_lengths.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped_basic_analyses.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped_lengths.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped_basic_analyses.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped_lengths.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped_basic_analyses.txt', sample=samples.index),
  		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped_lengths.txt', sample=samples.index),
#	R_script_files_from <plotting.R> script
		'Seq_Data_and_Output/Mapping_Distribution.txt',
#	organising outputfiles
		'Seq_Data_and_Output/Plots/Mapping_Distribution.txt'
##################################################################################################################################
########################             Adaptor_Trimming          ###################################################################
##################################################################################################################################
rule Adaptor_Trimming:
	input: 
		linker3 = 'refs/linker3.txt',
		samples = 'Seq_Data_and_Output/{sample}.fastq'
	output:
		trimmed = 'Seq_Data_and_Output/{sample}_Trimmed.txt',
		untrimmed = 'Seq_Data_and_Output/{sample}_Untrimmed.txt',
		trim_output = 'Seq_Data_and_Output/{sample}_Trimming_Output.txt'
	priority: 40
	shell:
		'tools/btrim32-static -p {input.linker3} -t {input.samples} -o {output.trimmed} -K {output.untrimmed} -3 -l 15 > {output.trim_output}'

##################################################################################################################################
########################             Qualtiy_Control           ###################################################################
##################################################################################################################################
rule Quality_Controll:
	input:
		samples = expand('Seq_Data_and_Output/{sample}_Trimmed.txt', sample=samples.index)
	output:
		fastqc = expand('Seq_Data_and_Output/FastQC/{sample}_Trimmed_fastqc.html', sample=samples.index)
	priority: 40
	shell:
		'fastqc {input.samples} -t 8 --noextract -o Seq_Data_and_Output/FastQC'
##################################################################################################################################
########################             rRNA_mapping              ###################################################################
##################################################################################################################################
rule rRNA_mapping:
	input: 
		samples = 'Seq_Data_and_Output/{sample}_Trimmed.txt',
	output:
		rRNA_unmapped = 'Seq_Data_and_Output/{sample}_rRNA_unmapped.txt',
		rRNA_readcount = 'Seq_Data_and_Output/{sample}_rRNA_readcount.txt'
	conda:
		'envs/bowtie.yaml'
	priority: 40
	shell:
		'bowtie -p 4 --best --norc -v 1 --chunkmbs 1024 -t tools/bowtie_references/rRNA {input.samples} --un {output.rRNA_unmapped}\
		| awk -f src/rRNA_readcount.awk | awk -f src/rRNA_readcount2.awk > {output.rRNA_readcount}'
##################################################################################################################################
########################             tRNA_mapping             ###################################################################
##################################################################################################################################
rule tRNA_mapping:
	input:
		samples = 'Seq_Data_and_Output/{sample}_rRNA_unmapped.txt'
	output:
		r_tRNA_unmapped ='Seq_Data_and_Output/{sample}_r_tRNA_unmapped.txt',
		tRNA_readcount = 'Seq_Data_and_Output/{sample}_tRNA_readcount.txt'
	conda:
		'envs/bowtie.yaml'
	priority: 40
	shell:
		'bowtie -p 4 --best --norc -v 1 --chunkmbs 1024 -t tools/bowtie_references/tRNA {input.samples} --un {output.r_tRNA_unmapped}\
		| awk -f src/rRNA_readcount.awk | awk -f src/rRNA_readcount2.awk > {output.tRNA_readcount}'
##################################################################################################################################
########################             isomir_mapping           ###################################################################
##################################################################################################################################
rule isomir_fasta_generation:
	input:
		hairpin='{}.fa'.format(hairpin),
		mature_db='{mature}.fa'
	output:
		'{mature}_isomir.fa'
	conda:
		'envs/iso_ref.yaml'
	priority: 40
	shell:
		'python src/create_isomir_ref.py -H {input.hairpin} -M {input.mature_db} -3T refs/3prime_template.json -5T refs/5prime_template.json -OUT {output} -MOUT multireads.txt'

rule isomir_fasta_splitting:
	input: temp('{mature}_isomir.fa')
	output: expand('{{mature}}_isomir_{length}nt.txt', length=lengths)
	priority: 40
	shell:
		'awk -f  src/fasta_file_splitting.awk {input}'


def get_prefix(wildcards):
	return('refs/bowtie_index_{}nt/{}_isomir'.format(wildcards.length,mature))

rule isomir_alignment_index:
	input:
		file = '{mature}_isomir_{length}nt.txt',
	output:
		reference1='refs/bowtie_index_{length}nt/{mature}_isomir.1.ebwt',
		reference2='refs/bowtie_index_{length}nt/{mature}_isomir.2.ebwt',
		reference3='refs/bowtie_index_{length}nt/{mature}_isomir.3.ebwt',
		reference4='refs/bowtie_index_{length}nt/{mature}_isomir.4.ebwt',
		referencerev1='refs/bowtie_index_{length}nt/{mature}_isomir.rev.1.ebwt',
		referencerev2='refs/bowtie_index_{length}nt/{mature}_isomir.rev.2.ebwt'
	params:
		prefix = lambda wildcards: get_prefix(wildcards)
	conda:
		'envs/bowtie.yaml'
	priority: 40
	shell:
		'bowtie-build -q {input} {params.prefix}'

rule fastq_splitting:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_unmapped.txt'
	output:
		temp(expand('Seq_Data_and_Output/{{sample}}_r_tRNA_unmapped_{length}nt.txt', length=lengths))
	priority: 40
	shell:
		'awk -f  src/fastq_file_splitting.awk {input.samples}'

rule bowtie_alignment:
	input: 
		file='Seq_Data_and_Output/{sample}_r_tRNA_unmapped_{length}nt.txt',
		reference1='refs/bowtie_index_{{length}}nt/{}_isomir.1.ebwt'.format(mature)
	output: 
		sam=temp('{sample}_isomir_{length}nt_alignment.sam'),
		unmapped=temp('{sample}_isomir_{length}nt_unmapped.fastq')
	params:
		prefix = lambda wildcards: get_prefix(wildcards)
	threads: 4
	conda:
		'envs/bowtie.yaml'
	priority: 40
	shell:
		'bowtie --quiet -p 4 --best --sam --norc -v 3 --chunkmbs 1024 -t {params.prefix}\
		{input.file} --un {output.unmapped} > {output.sam}'

rule unmapped_fastq_merging:
	input: expand('{{sample}}_isomir_{length}nt_unmapped.fastq', length=lengths)
	output: 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped.txt' 
	priority: 40
	shell:
		'cat {input} > {output}'

rule sam_processing_length_merging:
	input: expand('{{sample}}_isomir_{length}nt_alignment.sam',length=lengths)
	output: 'Seq_Data_and_Output/{sample}_isomir_reads.txt'
	priority: 40
	shell:
		'awk -f src/sam_processing_length_merging.awk {input} > {output}'

rule isomir_read_counting:
	input: 'Seq_Data_and_Output/{sample}_isomir_reads.txt'
	output: 'Seq_Data_and_Output/{sample}_isomir_readcount.txt'
	priority: 40
	shell:
		'awk -f src/isomir_read_counting.awk {input} > {output}'

rule count_merging:
	input: expand('Seq_Data_and_Output/{sample}_isomir_readcount.txt',sample=samples.index)
	output: 'Seq_Data_and_Output/isomir_readcount_all.txt'
	params: samples = lambda wildcards: samples.index
	conda:
		'envs/merge.yaml'
	priority: 40
	script:
		'src/count_merging.R'

#################################################################################################################################
########################             snRNA_mapping           ###################################################################
##################################################################################################################################
rule snRNA_mapping:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped.txt'
	output:
		r_tRNA_isomir_snRNA_unmapped ='Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped.txt',
		snRNA_readcount = 'Seq_Data_and_Output/{sample}_snRNA_readcount.txt'
	conda:
		'envs/bowtie.yaml'
	priority: 40
	shell:
		'bowtie -p 4 --best --norc -v 1 --chunkmbs 1024 -t tools/bowtie_references/snRNA {input.samples} --un {output.r_tRNA_isomir_snRNA_unmapped}\
		| awk -f src/rRNA_readcount.awk | awk -f src/rRNA_readcount2.awk > {output.snRNA_readcount}'
#################################################################################################################################
########################             snoRNA_mapping           ###################################################################
##################################################################################################################################
rule snoRNA_mapping:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped.txt'
	output:
		r_tRNA_isomir_sn_snoRNA_unmapped ='Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped.txt',
		snoRNA_readcount = 'Seq_Data_and_Output/{sample}_snoRNA_readcount.txt'
	conda:
		'envs/bowtie.yaml'
	priority: 40
	shell:
		'bowtie -p 4 --best --norc -v 1 --chunkmbs 1024 -t tools/bowtie_references/snoRNA {input.samples} --un {output.r_tRNA_isomir_sn_snoRNA_unmapped}\
		| awk -f src/rRNA_readcount.awk | awk -f src/rRNA_readcount2.awk > {output.snoRNA_readcount}'
#################################################################################################################################
########################            PLOTTING Info       ############################################################
##################################################################################################################################
rule trimming_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_Trimming_Output.txt'
	output:
		short_noadapter = 'Seq_Data_and_Output/{sample}_short_noadapter_1.txt'
	priority: 40
	shell:
		'grep -E "Short|No 3" {input.samples} > {output.short_noadapter}'

rule trimming_info_2:
	input:
		samples = 'Seq_Data_and_Output/{sample}_short_noadapter_1.txt'
	output:
		short_noadapter = 'Seq_Data_and_Output/{sample}_short_noadapter.txt'
	priority: 40
	shell:
		"sed '2~2 s/^.*:/No_Adapter:/' {input.samples} > {output.short_noadapter}"

rule length_and_seq_ab_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_Trimmed.txt'
	output:
		trimmedreads_basic_analyses = 'Seq_Data_and_Output/{sample}_trimmedreads_basic_analyses.txt',
		trimmedreads_lengths = 'Seq_Data_and_Output/{sample}_trimmedreads_lengths.txt'
	priority: 40
	shell:
		'perl src/TBr2_basic-analyses.pl -i {input.samples} -o {output.trimmedreads_basic_analyses} && awk -f src/trimmedreads_length.awk {output.trimmedreads_basic_analyses} > {output.trimmedreads_lengths}'


rule untrimmed_length_info:
	input:
		linker3 = 'refs/linker3.txt',
		samples = 'Seq_Data_and_Output/{sample}_Untrimmed.txt'
	output:
		short = 'Seq_Data_and_Output/{sample}_Short.txt',
		tmp   = 'Seq_Data_and_Output/{sample}_tmp.txt',
		short_lengths = 'Seq_Data_and_Output/{sample}_short_lengths.txt'
	priority: 40
	run:
		shell('tools/btrim32-static -p {input.linker3} -t {input.samples} -o {output.short} -3 -l 0 | src/cut_short.sh > {output.tmp};\
			cat {output.short} | src/short_length.sh | join {output.tmp} - > {output.short_lengths};')


rule rRNA_length_and_seq_ab_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_rRNA_unmapped.txt'
	output:
		rRNA_unmapped_basic_analyses = 'Seq_Data_and_Output/{sample}_rRNA_unmapped_basic_analyses.txt',
		rRNA_unmapped_lengths = 'Seq_Data_and_Output/{sample}_rRNA_unmapped_lengths.txt'
	priority: 40
	shell:
		'perl src/TBr2_basic-analyses.pl -i {input.samples} -o {output.rRNA_unmapped_basic_analyses} && awk -f src/rRNA_unmapped_basic_analyses.awk {output.rRNA_unmapped_basic_analyses} > {output.rRNA_unmapped_lengths}'


rule r_tRNA_length_and_seq_ab_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_unmapped.txt'
	output:
		r_tRNA_unmapped_basic_analyses = 'Seq_Data_and_Output/{sample}_r_tRNA_unmapped_basic_analyses.txt',
		r_tRNA_unmapped_lengths = 'Seq_Data_and_Output/{sample}_r_tRNA_unmapped_lengths.txt'
	priority: 40
	shell:
		'perl src/TBr2_basic-analyses.pl -i {input.samples} -o {output.r_tRNA_unmapped_basic_analyses} && awk -f src/rRNA_unmapped_basic_analyses.awk {output.r_tRNA_unmapped_basic_analyses} > {output.r_tRNA_unmapped_lengths}'

rule r_tRNA_isomir_length_and_seq_ab_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped.txt'
	output:
		r_tRNA_isomir_unmapped_basic_analyses = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped_basic_analyses.txt',
		r_tRNA_isomir_unmapped_lengths = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_unmapped_lengths.txt'
	priority: 40
	shell:
		'perl src/TBr2_basic-analyses.pl -i {input.samples} -o {output.r_tRNA_isomir_unmapped_basic_analyses} && awk -f src/rRNA_unmapped_basic_analyses.awk {output.r_tRNA_isomir_unmapped_basic_analyses} > {output.r_tRNA_isomir_unmapped_lengths}'

rule r_tRNA_isomir_snRNA_length_and_seq_ab_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped.txt'
	output:
		r_tRNA_isomir_snRNA_unmapped_basic_analyses = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped_basic_analyses.txt',
		r_tRNA_isomir_snRNA_unmapped_lengths = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_snRNA_unmapped_lengths.txt'
	priority: 40
	shell:
		'perl src/TBr2_basic-analyses.pl -i {input.samples} -o {output.r_tRNA_isomir_snRNA_unmapped_basic_analyses} && awk -f src/rRNA_unmapped_basic_analyses.awk {output.r_tRNA_isomir_snRNA_unmapped_basic_analyses} > {output.r_tRNA_isomir_snRNA_unmapped_lengths}'

rule r_tRNA_isomir_sn_snoRNA_length_and_seq_ab_info:
	input:
		samples = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped.txt'
	output:
		r_tRNA_isomir_sn_snoRNA_unmapped_basic_analyses = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped_basic_analyses.txt',
		r_tRNA_isomir_sn_snoRNA_unmapped_lengths = 'Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped_lengths.txt'
	priority: 40
	shell:
		'perl src/TBr2_basic-analyses.pl -i {input.samples} -o {output.r_tRNA_isomir_sn_snoRNA_unmapped_basic_analyses} && awk -f src/rRNA_unmapped_basic_analyses.awk {output.r_tRNA_isomir_sn_snoRNA_unmapped_basic_analyses} > {output.r_tRNA_isomir_sn_snoRNA_unmapped_lengths}'

#################################################################################################################################
########################            Run R script for plots      ############################################################
##################################################################################################################################


rule run_R_plotting_script:
	input:
		expand('Seq_Data_and_Output/{sample}_r_tRNA_isomir_sn_snoRNA_unmapped_basic_analyses.txt', sample=samples.index)
	output:
		'Seq_Data_and_Output/Mapping_Distribution.txt'
	priority: 1
	conda:
		'envs/r_plotting.yaml'
	script:
		'src/plotting.R'

##################################################################################################################################
#########################            Organizing files into folders     ############################################################
###################################################################################################################################
rule organising_files:
	input: 
		inp = 'Seq_Data_and_Output/Mapping_Distribution.txt'
	output:
		out = 'Seq_Data_and_Output/Plots/Mapping_Distribution.txt' 
	priority: 0
	shell:
		'mv Seq_Data_and_Output/*Distribution* Seq_Data_and_Output/Plots;'
		'mv Seq_Data_and_Output/*Overview* Seq_Data_and_Output/Plots/Additional_Files;'
		'mv Seq_Data_and_Output/*Readcount.txt Seq_Data_and_Output/Readcounts;'
		'mv Seq_Data_and_Output/*rRNA_readcount.txt Seq_Data_and_Output/Readcounts/rRNA;'
		'mv Seq_Data_and_Output/*snRNA_readcount.txt Seq_Data_and_Output/Readcounts/snRNA;'
		'mv Seq_Data_and_Output/*snoRNA_readcount.txt Seq_Data_and_Output/Readcounts/snoRNA;'
		'mv Seq_Data_and_Output/*tRNA_readcount.txt Seq_Data_and_Output/Readcounts/tRNA;'
		'mv Seq_Data_and_Output/*isomir_readcount.txt Seq_Data_and_Output/Readcounts/isomir;'
		'mv Seq_Data_and_Output/*isomir_unmapped.txt Seq_Data_and_Output/Unmapped;'
		'mv Seq_Data_and_Output/*Trimmed.txt Seq_Data_and_Output/Trimming;'
		'mv Seq_Data_and_Output/*Trimming_Output.txt Seq_Data_and_Output/Trimming/Trimming_Output;'
		'mv Seq_Data_and_Output/*Untrimmed.txt Seq_Data_and_Output/Trimming/Untrimmed;'
		'rm Seq_Data_and_Output/*txt;'