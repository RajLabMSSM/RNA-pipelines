# RAPiD-nf

RAPiD is an RNASeq pipeline used and maintained by Ying-Chih Wang and Hardik Shah in the Department of Genetics and Genomic Sciencese, Icahn Institute for Genomics and Multiscale Biology at Mount Sinai [1]. While RAPiD is built on Mount Sinai's APOLLO framework [2], RAPiD-nf inherits the same pipeline structure and is ported to Nextflow pipeline engine [3].

Authors
-------
 - Ying-Chih Wang
 - Hardik Shah

Collaborators
-------------
 - Kaur Alasoo
 - Gabriel Hoffman
 - Veera Manikandan
 - Panos Roussos
 - Solly Sieberts
 - Laura Sloofman
 - Georgios Voloudakis
 

Reference
---------
[1] Y.-C. Wang et al., "RAPiD: An Agile and Dependable RNA-Seq Framework," in The 65th Annual Meeting of The American Society of Human Genetics, Baltimore, MD, 2015.

[2] R. Z. Castellanos et al., "Apollo: A Production Tested, Vertically Integrated, Operations Enhanced, Science Aware Framework for Launching Large Cohorts of Genomics Pipelines," in The 65th Annual Meeting of The American Society of Human Genetics, Baltimore, MD, 2015.

[3] Nextflow Data-driven computational pipelines, https://www.nextflow.io/

License
-------
TBD

# Raj Lab stuff - added by Jack

Running on example data
-----------------------

Clone the RAPiD-nf repo

```
git clone git@github.com:genomely/RAPiD-nf.git
cd RAPiD-nf
```

RAPiD works on FASTQ files as input. Assume you have a list of FASTQ files with the format SAMPLE1_R1.fastq.gz, SAMPL1_R2.fastq.gz.

First symlink all FASTQ files to a directory called fastq/

```
for i in $(cat my_files.txt); do
    ln -s $i fastq/
done
```

Then run `create_fastq_key.R` - this creates a sample key for you with three columns - sample, f1, f2 (read 1 and read 2).


Then run `create_rapid_structure.py` - this creates the directory structure RAPiD-nf expects.
 
### Running on Chimera

First edit config/lsf.config to change the `acc_apollo` code to your account code, eg acc_ad-omics.

Then start a screen session:

```
screen -S RAPiD
```

You can run RAPiD like so:

```
bsub -I -n 2 -W 24:00 -q premium -P acc_als-omics -R rusage[mem=3750] -R span[hosts=1] "/sc/hydra/projects/PBG/nextflow/bin/nextflow run RAPiD.nf --run `pwd` --genome GRCh38.Gencode.v30 --stranded reverse -profile chimera --qc --fastqc --leafcutter --featureCounts --rsem --trimAdapter --rawPath Raw/Illumina -resume"
```

You can also submit to express instead of premium. Edit config/lsf.config to swap out premium for express.

``` bsub -I -n 2 -W 12:00 -q express -P acc_als-omics -R rusage[mem=3750] -R span[hosts=1] "/sc/hydra/projects/PBG/nextflow/bin/nextflow run RAPiD.nf --run `pwd` --genome GRCh38.Gencode.v30 --stranded reverse -profile chimera --qc --fastqc --leafcutter --featureCounts --rsem --trimAdapter --rawPath Raw/Illumina -resume"
```


Specifying Parameters
---------------------

Boolean values for which steps to run, as well as specifications about the RNA-seq samples, are stored in `experiment.config`. Setting `-profile experiment` in the nextflow command tells Nextflow to look inside this file, as well as `config/chimera.config` and `lsf.config`, which gives parameters for use on Chimera and parallel job submission using LSF. The use of these config files are set by `nextflow.config` and can be altered by the user.


Added by Katia Lopes: 

Checking strand      
---------------------

To read: http://onetipperday.sterding.com/2012/07/how-to-tell-which-library-type-to-use.html  

1 - Check the stranded with blat from UCSC 
Input: Part of sequence from R1 
Part of sequence from R2 
Look at the ID for the strand column and check if the same read is aligned to a different strand 

Taking a real example, first getting some reads (in fasta format) from the paired-end sequencing fastq file using command like:

cd /sc/arion/projects/als-omics/microglia_stimulated/raw_data/30-410768147/

zcat fastq/Sample_167_END1.fastq.gz | sed 's/@//g;s/ /_/g' | awk '{if(NR%4==1)print ">"$0;if(NR%4==2) print $0;}' | head
zcat fastq/Sample_167_END2.fastq.gz | sed 's/@//g;s/ /_/g' | awk '{if(NR%4==1)print ">"$0;if(NR%4==2) print $0;}' | head

https://genome.ucsc.edu/cgi-bin/hgBlat 

2 -  check the strandness of the reads with RSeQC: 
Run RAPID for few samples because you'll need a bam file for it. 
Run the script: infer_experiment.py (already installed on Chimera)

Example for test: 
infer_experiment.py -r ~/ENSEMBL_ANNOTATION/hg38_Gencode_V30.bed -i Sample_511/Processed/RAPiD/bams/Sample_511.bam


Checking adaptors
---------------------
Options: 

NexteraPE-PE 
TruSeq2-PE 
TruSeq2-SE 
TruSeq3-PE-2 
TruSeq3-PE 
TruSeq3-SE 

SE = Single End
PE = Paired End 

Run RAPID for few samples and check the percentage of adaptors. Adaptor from Microglia Stimulated samples: NexteraPE-PE

Command line
---------------------
bsub -I -n 2 -W 144:00 -q long -P acc_ad-omics -R rusage[mem=3750] -R span[hosts=1] "/sc/arion/projects/H_PBG/nextflow/bin/nextflow run /sc/arion/projects/als-omics/microglia_stimulated/RAPID_results/RAPiD-nf_1/RAPiD.nf --run `pwd` --genome GRCh38.Gencode.v30 --stranded none -profile chimera --qc --fastqc --leafcutter --featureCounts --rsem --kallisto --salmon --trimAdapter NexteraPE-PE --rawPath Raw/Illumina -resume"

