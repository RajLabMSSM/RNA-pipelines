module load sratoolkit
ml pigz
accessions=$1

mkdir -p fastq
for i in $(cat $accessions ) ; do 
	echo $i
    if [[ ! -f fastq/${i}_1.fastq && ! -f fastq/${i}_1.fastq.gz ]]; then
        prefetch $i --max-size u
        fasterq-dump --progress --include-technical --split-files -t tmp/ -O fastq/ -e 4 $i
    fi
    if [ ! -f fastq/${i}_1.fastq.gz ]; then
        pigz -p 4 fastq/${i}_1.fastq
        pigz -p 4 fastq/${i}_2.fastq
        pigz -p 4 fastq/${i}_3.fastq
        rm -r ${i}/
        rm -r tmp/
    fi
done
