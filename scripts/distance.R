# calculate distance between alignments using Blackburne and Whelan 2012
# doi: https://doi.org/10.1093/bioinformatics/btr701
# binary `metal` calculates the distance available on GitHub by the authors

distance = function(filename, models, metric, brlen) {
    # choose between d_seq and d_pos
    if(metric == "dseq") {
        flag = paste("-s")
    } else if(metric == "dpos") {
        flag = paste("-p")
    } else {
        stop(paste("Metric", metric, "not supported."))
    }
    
    ref = paste0("results/", brlen, "/ref_alignments/", filename)
    distances = c()
    # calculate distance between reference for each method
    for(model in models) {
        aln = paste0("results/", brlen, "/aln/", model, "/", filename)
        # call metal to calculate distance
        d = suppressWarnings(system(paste("bin/metal", flag, ref, aln), intern = TRUE, ignore.stderr = TRUE))
        # parse metal output ("fraction = double" --> double)
        if(length(d) <= 0) {
            d = NA
        } else {
            d = as.double(gsub("^.*?=\ ", "", d))
        }
        distances = append(distances, d)
    }
    
    return(distances)
}
