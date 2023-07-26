# filter sequences by length - discard long sequences

N=$1
SIZE=$2
OUTPUT=$3

rm -f ${OUTPUT}

for file in $(ls data/*.fasta | head -n${N} )
do
	s=$(stat --printf="%s" ${file})
	if [ ${s} -le ${SIZE} ]
	then
		echo $(basename ${file}) >> ${OUTPUT}
	fi
done

if [ ! -f ${OUTPUT} ]
then
	echo "No files smaller than ${SIZE} found."
fi
