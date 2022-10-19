bsub -I -n 2 -W 144:00 -q long -P acc_als-omics -R rusage[mem=3750] -R span[hosts=1] "/sc/arion/projects/H_PBG/nextflow/bin/nextflow run RAPiD.nf --run `pwd` --genome GRCh38.Gencode.v30 --stranded none -profile TwoPassSTAR --qc --fastqc --leafcutter --featureCounts --rsem --kallisto --salmon --rawPath Raw/Illumina --trimAdapter TruSeq3-PE"
