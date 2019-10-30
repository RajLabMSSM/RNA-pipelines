#!/bin/sh
# properties = {"type": "single", "rule": "multiQC", "local": false, "input": [], "output": ["/sc/orga/projects/als-omics/NYGC_ALS/data/oct_2019_multiqc_report.html"], "wildcards": {}, "params": {"out": "/sc/orga/projects/als-omics/NYGC_ALS/data/oct_2019_"}, "log": [], "threads": 1, "resources": {}, "jobid": 2, "cluster": {"queue": "premium", "cores": 1, "mem": 32000, "time": "180", "name": "$(basename $(pwd)):multiQC:", "output": "logs/multiQC:.stdout", "error": "logs/multiQC:.stderr"}}
cd /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples && \
/sc/orga/work/humphj04/conda/envs/leafcutter-pipeline/bin/python3.6 \
-m snakemake /sc/orga/projects/als-omics/NYGC_ALS/data/oct_2019_multiqc_report.html --snakefile /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /sc/orga/projects/als-omics/NYGC_ALS/pipelines/collate_samples/.snakemake/tmp.i_xcya53 --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules multiQC --nocolor --notemp --no-hooks --nolock \
--mode 2 

