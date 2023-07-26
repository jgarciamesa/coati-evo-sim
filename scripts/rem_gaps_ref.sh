# create a file with the name of the pairwise aligned ENSEMBL sequences
#  that DON'T contain gaps

for file in $(ls results/$1/ref_alignments/*)
do
	cat ${file} | tr -d '-' > results/$1/no_gaps_ref/$(basename ${file})
done

