#!/bin/sh
# properties = {"type": "single", "rule": "genes_pca", "local": false, "input": ["/sc/orga/projects/als-omics/RAPiD/gene_matrix.RData", "/sc/orga/projects/als-omics/NYGC_ALS/tables/rna_support.tsv"], "output": ["/sc/orga/projects/als-omics/NYGC_ALS/data/genes_pca.RData"], "wildcards": {}, "params": {"script": "/sc/orga/projects/als-omics/NYGC_ALS/scripts/genes_pca.R"}, "log": [], "threads": 1, "resources": {}, "jobid": 2, "cluster": {"queue": "premium", "cores": 1, "mem": 7750, "time": 60, "name": "$(basename $(pwd)):genes_pca:", "output": "logs/genes_pca:.stdout", "error": "logs/genes_pca:.stderr"}}
cd /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples && \
/sc/orga/work/humphj04/conda/envs/isoseq-pipeline/bin/python3.6 \
-m snakemake /sc/orga/projects/als-omics/NYGC_ALS/data/genes_pca.RData --snakefile /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples/.snakemake/tmp.nmpkxzer /sc/orga/projects/als-omics/RAPiD/gene_matrix.RData /sc/orga/projects/als-omics/NYGC_ALS/tables/rna_support.tsv --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules genes_pca --nocolor --notemp --no-hooks --nolock \
--mode 2 

