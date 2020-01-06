#ml bioawk
bioawk=/hpc/packages/minerva-common/bioawk/1.0/bin/bioawk
if [ ! -e $1 ];then
	echo "file doesn't exist"
fi
zless $1 | $bioawk -c fastx '{print $name}' | wc -l 
