#!/bin/sh
# properties = {"type": "single", "rule": "getMatrices", "local": false, "input": [], "output": ["/sc/orga/projects/als-omics/RAPiD/gene_matrix.RData", "/sc/orga/projects/als-omics/RAPiD/tx_matrix.RData"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 1, "cluster": {"queue": "premium", "cores": 1, "mem": 3750, "time": "180", "name": "$(basename $(pwd)):getMatrices:", "output": "logs/getMatrices:.stdout", "error": "logs/getMatrices:.stderr"}}
cd /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples && \
/sc/orga/work/humphj04/conda/envs/isoseq-pipeline/bin/python3.6 \
-m snakemake /sc/orga/projects/als-omics/RAPiD/gene_matrix.RData --snakefile /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples/.snakemake/tmp.qu2dx8w_ --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules getMatrices --nocolor --notemp --no-hooks --nolock \
--mode 2 

