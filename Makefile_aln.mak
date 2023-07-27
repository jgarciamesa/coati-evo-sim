################################################################################
# aligner recipes for coati: tri-mg, mar-mg, tri-ecm, mar-ecm				   #
################################################################################

results/$(T)/aln/tri-mg/%: results/$(T)/no_gaps_ref/%
	@echo -ne "coati tri-mg align $*\t\t\r            "
	@./bin/coati-alignpair $< -m tri-mg -o $@ -t $(T)

results/$(T)/aln/mar-mg-sum/%: results/$(T)/no_gaps_ref/%
	@echo -ne "mcoati mar-mg-sum align $*\t\t\r       "
	@./bin/coati-alignpair $< -m mar-mg -o $@ -t $(T) --marginal-sub SUM

results/$(T)/aln/mar-mg-max/%: results/$(T)/no_gaps_ref/%
	@echo -ne "mcoati mar-mg-max align $*\t\t\r       "
	@./bin/coati-alignpair $< -m mar-mg -o $@ -t $(T) --marginal-sub MAX

clean_aln:
	rm -f results/$(T)/aln/tri-mg/*.fasta
	rm -f results/$(T)/aln/mar-mg-sum/*.fasta
	rm -f results/$(T)/aln/mar-mg-max/*.fasta

