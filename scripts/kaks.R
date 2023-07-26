suppressMessages(library(seqinr))
library(stringr)

#################### code modified from Reed A Cartwright ######################
ka = function(a) { if(class(a) != "alignment") {return(NA)}
	ka = kaks(a)$ka[1]
	ka
}

ks = function(a) {
	if(class(a) != "alignment") {return(NA)}
	ks = kaks(a)$ks[1]
	ks
}

omega = function(a) { #dN/dS
    if(class(a) != "alignment") {return(NA)}
	k = kaks(a)
	ka = k$ka[1]
	ks = k$ks[1]
	if(ka < 0 || ks < 0) {
		return(0)
	} else if(ka == 0 && ks == 0) {
		return(0)
	} else if(ks == 0) {
		return(10)
	} else {
		return(ka/ks)
	}
}

ps_accuracy = function(ref,w) {
	tp = length(intersect(which(w>1),which(ref>1)))
	fp = length(intersect(which(w>1),which(ref<1)))
	fn = length(intersect(which(w<1),which(ref>1)))
	if(all(c(tp,fp,fn) == 0)) {
	    return(NA)
	}
	return(2*tp/(2*tp+fp+fn))
}

ns_accuracy = function(ref,w) {
	tp = length(intersect(which(w<1),which(ref<1)))
	fp = length(intersect(which(w<1),which(ref>1)))
	fn = length(intersect(which(w>1),which(ref<1)))
	if(all(c(tp,fp,fn) == 0)) {
	    return(NA)
	}
	return(2*tp/(2*tp+fp+fn))
}

process_aln = function(d,fasta) {
	ret = list()

	for(f in fasta) {
		aln = try(read.alignment(str_c(d,"/",f),"fasta"),silent=TRUE)
		if(class(aln) == "try-error") {
			ret[[f]] = list()
		} else {
			seq = aln$seq %>% str_split(pattern="")
			ung = (seq[[1]] %in% c("a","c","g","t","n"))
			seq = lapply(seq,subset,ung) # subset = "["
			aln$seq = lapply(seq,str_c,collapse="")
			ret[[f]] = aln
		}
	}
	ret
}

################################################################################

k_metric = function(ARGS) {
    results = data.frame("ka RMSE" = double(),
                         "ks RMSE" = double(),
                         "F1-score pos selection" = double(),
                         "F1-score neg selection" = double())
	for(i in ARGS){
		# read alignments
		model = list.files(paste0("aln/ref/",i), pattern = "*.fasta")
		ref = list.files(paste0("data/ref_alignments"),pattern = "*.fasta")

		if(length(ref) != length(model)) {
			warning(paste("Number of sequences between ref and", i, "differs by", abs(length(ref)-length(model))))
		}

		# remove insertions w.r.t. model organism (human)
		ref = process_aln(paste0("data/ref_alignments"),model)
		model = process_aln(paste0("aln/ref/",i),model)

		# calculate ka
		ref_ka = sapply(ref,ka)
		model_ka = sapply(model,ka)

		# calculate ks
		ref_ks = sapply(ref,ks)
		model_ks = sapply(model,ks)

		# calculate omega
		ref_omega = sapply(ref,omega)
		model_omega = sapply(model,omega)

		# calculate root mean-square error of ka
		model_rmse_ka = sqrt(sum((ref_ka-model_ka)^2)/length(ref_ka))

		# calculate root mean-square error of ks
		model_rmse_ks = sqrt(sum((ref_ks-model_ks)^2)/length(ref_ks))

		# calculate accuracy of positive selection
		model_pos = ps_accuracy(ref_omega,model_omega) 

		# calculate accuracy of negative selection
		model_neg = ns_accuracy(ref_omega,model_omega) 

		results[i,] = c(model_rmse_ka, model_rmse_ks, model_pos, model_neg)
	}
    return(results)
}
