# the following must be set by the user
BASE = '/home/users/hjp/intron-retention'
GEUVADIS_DATA = '/scratch/PI/pritch/geuvadis'

## software

IR_BASE = '/home/users/hjp/kma/inst/pre-process'
FASTQTL = 'module load gsl/1.16 && /home/users/hjp/software/FastQTL/bin/fastQTL.static'

# end customizations

# INDEX = BASE
INDEX = BASE + '/index'

HUMAN_BASE = INDEX + '/Homo_sapiens.GRCh37.75.dna.primary_assembly'
HUMAN_GENOME = INDEX + '/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
HUMAN_GTF = INDEX + '/Homo_sapiens.GRCh37.75.gtf'

HUMAN_BOWTIE = INDEX + '/transcripts_and_introns'
HUMAN_KALLISTO_INDEX = INDEX + '/transcripts_and_introns.kidx'

MERGED_GTF = INDEX + '/merged.gtf'

def source_r(base, fname):
    return 'OMP_NUM_THREADS=1 module load R && cd {0} && R CMD BATCH --vanilla --default-packages=methods,stats,utils {1}'.format(base, fname)
