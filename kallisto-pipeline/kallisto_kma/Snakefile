include: 'config.py'
from os import listdir
from os.path import basename

import re

## JH: in future hardcode paths to human GTF and Genome for hg38
#HUMAN_GTF=""
#HUMAN_GENOME=""

shell.prefix("ml kallisto/0.45.0")

rule all:
    input:
        HUMAN_GTF,
        HUMAN_GENOME,
        expand(HUMAN_BASE + '.{i}.ht2', i = range(1, 9)),


        'kma_out/introns.fa',
        INDEX + '/transcripts.fa',
        INDEX + '/transcripts_and_introns.fa',

        # XXX: this will break right here
        # need to run `src/annotate_gtf.R` after running the stuff above
        expand(OUTFOLDER + '/results/kma/{sample}_genome/pseudoalignments.bam',
            sample = GEUVADIS_FILES),

rule get_human_gtf:
    output:
        HUMAN_GTF
    params:
        base_name = 'Homo_sapiens.GRCh37.75.gtf'
    shell:
        'wget -O {HUMAN_GTF}.gz ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz'
        ' && '
        'gunzip {params.base_name}.gz'

rule get_human_genome:
    output:
        HUMAN_GENOME
    params:
        base_name = 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
    shell:
        'wget -O {HUMAN_GENOME}.gz ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'
        ' && '
        'gunzip {params.base_name}.gz'

rule build_introns:
    input:
        HUMAN_GTF,
        HUMAN_GENOME
    output:
        BASE + '/kma_out',
        BASE + '/kma_out/introns.fa'
    shell:
        'python2 {IR_BASE}/generate_introns.py'
        ' --extend 50'
        ' --gtf {HUMAN_GTF}'
        ' --out {output[0]}'
        ' --genome {HUMAN_GENOME}'
        ' > '
        'kma.log'

rule get_transcriptome:
    input:
        HUMAN_GTF,
        HUMAN_GENOME
    output:
        INDEX + '/transcripts.fa'
    shell:
        'gffread'
        ' {HUMAN_GTF}'
        ' -g {HUMAN_GENOME}'
        ' -w {output}'

rule merge_transcriptome:
    input:
        INDEX + '/transcripts.fa',
        BASE + '/kma_out/introns.fa'
    output:
        INDEX + '/transcripts_and_introns.fa'
    shell:
        'cat {input[0]} {input[1]} > {output}'

rule kallisto_index:
    benchmark:
        'benchmark/kallisto_index.json'
    input:
        BASE + '/transcripts_and_introns.fa'
    output:
        HUMAN_KALLISTO_INDEX
    shell:
        'kallisto index'
        ' -i {output}'
        ' {input}'

rule kallisto_quant:
    benchmark:
        'benchmark/kallisto/{sample}.json'
    input:
	#JH - make this more general
        OUTFOLDER + '/fastq/{sample}_1.fastq.gz',
        OUTFOLDER + '/fastq/{sample}_2.fastq.gz',
        HUMAN_KALLISTO_INDEX
    output:
        OUTFOLDER + '/results/kma/{sample}',
        OUTFOLDER + '/results/kma/{sample}/abundance.h5'
    threads: 5
    shell:
        'kallisto quant'
        ' -i {HUMAN_KALLISTO_INDEX}'
        ' -o {output[0]}'
        ' -t {threads}'
        ' {input[0]} {input[1]}'

rule kallisto_genome:
    benchmark:
        'benchmark/kallisto_genome/{sample}.json'
    input:
        OUTFOLDER + '/fastq/{sample}_1.fastq.gz',
        OUTFOLDER + '/fastq/{sample}_2.fastq.gz',
        HUMAN_KALLISTO_INDEX
    output:
        OUTFOLDER + '/results/kma/{sample}_genome',
        OUTFOLDER + '/results/kma/{sample}_genome/abundance.h5',
        OUTFOLDER + '/results/kma/{sample}_genome/pseudoalignments.bam'
    threads: 5
    shell:
        'kallisto quant'
        ' -i {HUMAN_KALLISTO_INDEX}'
        ' -o {output[0]}'
        ' -t {threads}'
        ' --genomebam'
        ' --gtf {MERGED_GTF}'
        ' {input[0]} {input[1]}'
