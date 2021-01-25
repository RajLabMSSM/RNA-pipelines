# SNAKEJOB - running snakemake pipelines on Minerva

Assumes you have a cluster.yaml file written for your pipeline. See example file.

Each rule must have resources listed in cluster.yaml.

# USAGE

Both scripts assume you have a conda environment that contains snakemake, called "snakemake".

Check whether you do by running `conda activate snakemake`

Make snakejob and snakejob_HPC executable and findable on your `$PATH`

```
chmod +x snakejob
chmod +x snakejob_HPC
ln -s $PWD/snakejob /usr/bin/
ln -s $PWD/snakejob_HPC /usr/bin/
```

## Running options

0. snakemake runs locally and runs jobs locally

```
snakemake -s Snakefile --configfile config.yaml
```

1. snakemake runs locally but submits jobs to cluster

```
./snakejob -s Snakefile -c config.yaml

Usage: ./snakejob -c config.yaml -s Snakefile
Options:
    -c config file (default: config.yaml)
    -s Snakefile (default: Snakefile)
    -n dry run mode
    -m mode (if your pipeline requires a mode variable using --config)
    -a which HPC account to run on (default: acc_als-omics)
```

2. snakemake runs as a job on the cluster and submits further jobs to the cluster

```
./snakejob_HPC -s Snakefile -c config.yaml
Options:
    -c config file (default: config.yaml)
    -s Snakefile (default: Snakefile)
    -a cluster account (default: acc_als-omics)
    -q queue to submit to (default: premium)
    -t time of submission (default: 24:00)
    -m mode (if your Snakefile takes a mode value with --config)
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
