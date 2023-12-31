################################################################################
# Modify values in this section                                                #
################################################################################

# n: number of sequences
n=16000

# len: max length of sequence
len=3000

# j: number of parallel jobs
j=10

# t: branch length used for simulations
t=0.4

################################################################################
# Pipeline                                                                     #
################################################################################

# Download gene ID list
tput setaf 11; echo "Download gene ids                "
tput setaf 15
time make download_geneid T=${t} N=${n}
tput setaf 11; echo "Gene ids downloaded              "
tput setaf 15

# Download sequences
tput setaf 11; echo "Download sequences               "
tput setaf 15
time make download_genes N=${n} -j4 T=${t}
tput setaf 11; echo "Sequences downloaded             "
tput setaf 15

# Simulate benchmark alignments
tput setaf 11; echo "Simulate true alignments         "
tput setaf 15
time make results/${t}/ref_alignments.csv -j${j} N=${n}
tput setaf 11; echo "True alignments simulated        "
tput setaf 15

# Remove gaps from benchmark alignments
tput setaf 11; echo "Remove gaps from true data set   "
tput setaf 15
time make no_gaps_reference T=${t} N=${n}
tput setaf 11; echo "Gaps from true data set removed  "
tput setaf 15

# Align
tput setaf 11; echo "Align sequences"
tput setaf 15
time make align T=${t} -j${j} -i N=${n}
tput setaf 11; echo "Done aligning sequences"
tput setaf 15

# Result summary statistics
tput setaf 11; echo "Compute result summary statistics"
tput setaf 15
time make results/${t}/results_summary.csv -j${j} -i N=${n}
make results/${t}/gap_stats/{freq,phase,pos}-reference.csv N=${n}
make results/${t}/gap_stats/{freq,phase,pos}-tri-mg.csv N=${n}
make results/${t}/gap_stats/{freq,phase,pos}-mar-mg-sum.csv N=${n}
make results/${t}/gap_stats/{freq,phase,pos}-mar-mg-max.csv N=${n}
time make results_table_dseq T=${t} N=${n}
tput setaf 11; echo "Summary statistics computed      "
tput setaf 15
