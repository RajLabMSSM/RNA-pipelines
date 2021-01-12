# SNAKEJOB - running snakemake pipelines on Minerva

Assumes you have a cluster.yaml file written for your pipeline. See example file.

Each rule must have resources listed in cluster.yaml.


## Running options

0. snakemake runs locally and runs jobs locally

```
snakemake -s Snakefile --configfile config.yaml
```

1. snakemake runs locally but submits jobs to cluster

```
./snakejob -s Snakefile -c config.yaml -m <mode> 
```

2. snakemake runs as a job on the cluster and submits further jobs to the cluster

```
./snakejob_HPC -s Snakefile -c config.yaml -m <mode>
```

<mode> is optional, this is used by the QTL pipeline to set a config parameter when running pipeline, rather than through the config file.

3. as before but for running multiple snakemake runs with different config files

```
sh run_all_QTLs_parallel.sh
``` 

This is currently set up for the QTL pipeline when running multiple datasets in parallel for maximum efficiency.


## Relationships


snakejob calls snakemake with the user's snakefile and config.yaml set using the flags, then uses the cluster.yaml file to submit each rule as a job on minerva.

snakejob_HPC is a wrapper for submitting snakejob as its own job on the cluster.
