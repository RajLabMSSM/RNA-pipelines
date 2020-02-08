# take in metadata and create the directory structure for RAPiD
# Jack Humphrey


# first argument to script should be a table with 2+ columns:
# sample: a list of sample names
# f1: full paths to the read 1 FASTQ, must have .fastq.gz at the end
# f2: full paths to read 2 FASTQ, must have .fastq.gz at the end

import pandas as pd
import os
import sys

# usage
# python create_rapid_structure.py samples.tsv 

# RAPiD-nf looks for FASTQ files in a <sample/Raw/Illumina/ folder 
# if paired-end, R1 files must have R1 in the file name and R2 must have R2. 
# SRA files don't do this so when making the symlink I append R1 or R2

metadata = sys.argv[1]


exists = os.path.isfile(metadata)
if not exists:
    print( metadata , "does not exist!" )
    exit()

print("reading", metadata)
meta = pd.read_csv(metadata, sep = '\t')

# check column names
mode = False
if 'sample' in list(meta) and 'f1' in list(meta):
    mode = "single"
if 'sample' in list(meta) and 'f1' in list(meta) and 'f2' in list(meta):
    mode = "paired"

if mode == False:
    print( metadata, 'does not contain correct columns \"sample\" and \"f1\"' )

print("samples are %s-end" % mode)

samples = meta['sample']
f1 = meta['f1']

# check fastqs exist
f1_exists = [os.path.isfile(f) for f in f1]

if not all(f1_exists):
    print( "missing files!" )
    exit()      

if mode == "paired":
    f2 = meta['f2']

# iterate through the samples
for i in range(len(samples)):
    print("sample:", samples[i])
    # make a directory for each sample
    # make <sample>/Raw/Illumina/
    sample_dir = samples[i] + "/Raw/Illumina/"
    os.makedirs(sample_dir, exist_ok = True)
    # symlink f1 (and f2 if exists)
    f1_source = os.path.abspath(f1[i])
    f1_target = sample_dir + os.path.basename(f1[i]).split('.fastq.gz')[0] + '_R1' + ".fastq.gz"
    if not os.path.exists(f1_target):
        print(f1_target)
        os.symlink(f1_source, f1_target )   
    if mode == 'paired':
        f2_source = os.path.abspath(f2[i])
        f2_target = sample_dir + os.path.basename(f2[i]).split('.fastq.gz')[0] + '_R2' + ".fastq.gz"
        if not os.path.exists(f2_target):
            print(f2_target)
            os.symlink(f2_source, f2_target )

