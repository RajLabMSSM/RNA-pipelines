# Adapter trimming with Cutadapt

RAPID-nf already has trimmomatic implemented, but for certain projects we need more flexibility with adapter trimming, such as weird fron end adapters used by the Revelo low-input RNA kit.

This pipeline uses cutdapt, reference here: https://cutadapt.readthedocs.io/en/stable/reference.html

The sample metadata table can be created using the RAPiD-tips script `create_fastq_key.R`
