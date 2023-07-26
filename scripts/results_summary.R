

source("scripts/distance.R")
source("scripts/kaks.R")

main = function(filename, brlen) {
    # read info about reference alns creation
    models = c("tri-mg", "mar-mg-sum", "mar-mg-max")
    reference_df = read.csv(paste0("results/", brlen, "/ref_alignments.csv"))
    len = length(models)

    ############################################################################
    # calculate dseq & dpos for inferred alignments
    dseqs = distance(filename, models, "dseq", brlen)
    dposs = distance(filename, models, "dpos", brlen)
    
    ############################################################################
    # omega for reference aln
    ref = process_aln(paste0("results/", brlen, "/ref_alignments/"), filename)
    ref_omega = sapply(ref,omega)
    
    # omega for inferred alns
    alns = sapply(models, function(x) {process_aln(paste0("results/", brlen, "/aln/", x, "/"), filename)})
    models_omega = sapply(alns, omega)
    
    # score alignments
    scores = c()
    for(model in models) {
        s = suppressWarnings(system(paste0("bin/coati-alignpair -s results/", brlen, "/aln/", model, "/", filename), intern = TRUE, ignore.stderr = TRUE))
        scores = c(scores, as.double(s))
    }

    if(length(scores) != len) {
        scores = rep(NA, len)
    }

    r_score = as.double(suppressWarnings(system(paste0("bin/coati-alignpair -s results/", brlen, "/ref_alignments/", filename), intern = TRUE, ignore.stderr = TRUE)))
    if(length(r_score) == 0) {
        r_score = NA
    }
                                                #, intern = TRUE, ignore.stderr = TRUE)))
    # group distance (dseq) and selection (omega) results
    distance_selection_df = data.frame(ensembl_id = rep(filename, len),
                            model = models,
                            dseq = dseqs,
                            dpos = dposs,
                            ref_omega = rep(ref_omega, len),
                            aln_omega = models_omega,
                            ref_score = rep(r_score, len),
                            score = scores)

    # merge results metrics with gap origin information
    results = merge(reference_df, distance_selection_df, by = "ensembl_id")
    write.table(results,
              file = paste0("results/", brlen, "/results_summary.csv"),
              append = TRUE,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              sep = ",")
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly = TRUE)
    main(ARGS[1], ARGS[-1])
}