################################################################################
# Global variables                                                             #
################################################################################
RSCRIPT = Rscript --vanilla
SHELL = /bin/bash

N ?= 16000
LEN ?= 3000
T ?= 0.4

$(shell mkdir -p data/)
$(shell mkdir -p results/$(T)/{aln,figures,gap_stats,no_gaps_ref,ref_alignments} 2> /dev/null)
$(shell mkdir -p results/$(T)/aln/{tri-mg,mar-mg-sum,mar-mg-max} 2> /dev/null)

################################################################################
# Download gene ID table from Ensembl (use to update files)                    #
################################################################################

.PHONY: download_geneid
download_geneid: | data/human_cds_geneId.tsv

data/human_cds_geneId.tsv: scripts/get_geneId.R
	@echo "Downloading" $* "gene IDs"
	@$(RSCRIPT) $< $* $@

################################################################################
# Download first N sequences from ENSEMBL                                      #
################################################################################
GENES = $(addsuffix .fasta,$(shell head -n${N} data/human_cds_geneId.tsv))
DOWNLOAD_GENES = $(addprefix data/,$(GENES))

.PHONY: download_genes
download_genes: $(DOWNLOAD_GENES)

data/%.fasta: | scripts/get_sequences.R
	@$(RSCRIPT) scripts/get_sequences.R data/human_cds_geneId.tsv $* $(@D)
	@echo -ne "Downloaded " $* "("$(shell ls data/*.fasta | wc -l 2> /dev/null)" out of ${N}).\r"

################################################################################
# Filter sequences by length												   #
################################################################################

results/$(T)/filtered.csv: $(DOWNLOAD_GENES) scripts/filter_seqs.sh
	@bash scripts/filter_seqs.sh $(N) $(LEN) $@

################################################################################
# Simulate alignments using coati's triplet model							   #
################################################################################
SIM = $(addprefix results/$(T)/ref_alignments/,$(shell cat results/$(T)/filtered.csv | head -n$(N) 2> /dev/null))

results/$(T)/ref_alignments.csv: $(SIM) | results/$(T)/filtered.csv
	@echo "Done creating reference alignments                   "
	@sed -i '1 i\ensembl_id,cigar,brlen,omega' $@

SIM_R = scripts/biological_simulation.R scripts/MG94.R scripts/write_fasta.R

results/$(T)/ref_alignments/%.fasta: $(SIM_R)
	@echo -ne "Creating reference alignment $*\r"
	@timeout 20s $(RSCRIPT) $< data/$*.fasta $@ 0.4 | cut -d '"' -f 2 >> results/$(T)/ref_alignments.csv

################################################################################
# Create reference alignments with no gaps for testing                         #
################################################################################

.PHONY: no_gaps_reference
no_gaps_reference:
	@bash scripts/rem_gaps_ref.sh $(T)

################################################################################
# Align simulated alignments												   #
################################################################################
ALN = $(shell ls results/$(T)/no_gaps_ref/ 2> /dev/null)

ALIGN = no_gaps_reference \
 				$(addprefix results/$(T)/aln/mar-mg-sum/,$(ALN)) \
 				$(addprefix results/$(T)/aln/mar-mg-max/,$(ALN)) \
 				$(addprefix results/$(T)/aln/tri-mg/,$(ALN))

.PHONY: align
align: $(ALIGN)

################################################################################
# Summary statistics                                                           #
################################################################################
ALL_ALNS = $(shell ls results/$(T)/ref_alignments/ 2> /dev/null)
RES = $(addprefix results/$(T)/, $(addsuffix .res, $(ALL_ALNS)))

results/$(T)/results_summary.csv: $(RES)
	@sed -i '1 i\ensembl_id,cigar,brlen,w,model,dseq,dpos,ref_omega,aln_omega,ref_score,score' $@
	@echo -ne "Built $@   \n"

results/$(T)/%.res: scripts/results_summary.R
	@echo -ne "summary stats $*\r"
	@$(RSCRIPT) $< $* $(T)

.PHONY: results_table
results_table: scripts/results_report.R | results/$(T)/results_summary.csv
	@$(RSCRIPT) $< results/$(T)/results_summary.csv

results/$(T)/gap_stats/freq-reference.csv: bin/sasi
	@$< gap frequency results/$(T)/ref_alignments/*.fasta -o $@

results/$(T)/gap_stats/phase-reference.csv: bin/sasi
	@$< gap phase results/$(T)/ref_alignments/*.fasta -o $@

results/$(T)/gap_stats/pos-reference.csv: bin/sasi
	@$< gap position results/$(T)/ref_alignments/*.fasta -o $@

results/$(T)/gap_stats/freq-%.csv: bin/sasi
	@$< gap frequency results/$(T)/aln/$*/*.fasta -o $@

results/$(T)/gap_stats/phase-%.csv: bin/sasi
	@$< gap phase results/$(T)/aln/$*/*.fasta -o $@

results/$(T)/gap_stats/pos-%.csv: bin/sasi
	@$< gap position results/$(T)/aln/$*/*.fasta -o $@

results/$(T)/figures/gap-%.pdf: scripts/plot_gap_stats.R
	@$(RSCRIPT) $< $* $(T)
################################################################################
# Clean pipeline results except gene id list and raw fasta downloads		   #
################################################################################

.PHONY: clean_pipeline

clean_pipeline:
	rm -f results/$(T)/filtered.csv
	rm -f results/$(T)/ref_alignments.csv
	rm -f results/$(T)/no_gaps_ref/*.fasta
	rm -f results/$(T)/ref_alignments/*.fasta
	rm -f results/$(T)/aln/{tri-mg,mar-mg-sum,mar-mg-max}/*.fasta
	rm -f results/$(T)/results_summary.csv

include Makefile_aln.mak
