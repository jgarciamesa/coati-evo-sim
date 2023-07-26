summary_table = function(results_file) {
    # read results_summary data
    results = read.csv(results_file, header = TRUE)
    
    models = c("tri-mg", "mar-mg-sum", "mar-mg-max")
    # extract dseq values for each aligner
    dseq = data.frame(row.names = "$d_{seq}")
    for(model in models) {
        dseq[model] = mean(results$dseq[results$model == model], na.rm = TRUE)
    }

    ############################################################################
    # perfect, best, and imperfect alignments
    source("scripts/number_alignments.R")
    alns = num_alns(results_file)
    
    # fix row and column names and add to stats table
    colnames(alns) = c("Perfect alignments", "Best alignments", "Imperfect alignments")
    stats_table = rbind(dseq, t(alns))
    
    ############################################################################
    # selection (ka ks)
    source("scripts/kaks.R")
    
    selection = data.frame("F1.pos" = double(),
                           "F1.neg" = double())
    
    # calculate F1 score for positive and negative selection for each aligner
    for(model in models) {
        model_rows = results[results$model == model, ]
        pos = ps_accuracy(model_rows$ref_omega, model_rows$aln_omega)
        neg = ns_accuracy(model_rows$ref_omega, model_rows$aln_omega)
        
        selection[nrow(selection) + 1, ] = c(pos, neg)
    }
    
    # fix row and column names and add to stats table
    colnames(selection) = c("F1-score pos selection", "F1-score neg selection")
    rownames(selection) = models
    stats_table = rbind(stats_table, t(selection))
    
    # display stats table
    stats_table
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly = TRUE)
    if(length(ARGS) != 1) {
        stop("Wrong number of arguments, please specify the results file.")
    }
    print(summary_table(ARGS[1]))
}