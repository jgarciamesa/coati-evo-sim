library(seqinr)
library(stringr)

source("scripts/MG94.R")

# Encode a sequence as a vector of integers where each codon is translated to
#  its position in the codon table (AAA -> 1, AAC -> 2, ... , TTT -> 61).
#  Stop codons are encoded as 0.
encode_seq = function(seqA) {
    stopifnot(nchar(seqA) %% 3 == 0)
    codons = get_codonstr(61)
    result = unlist(lapply(X = unlist(strsplit(gsub("(.{3})", "\\1 ", seqA), split = " ")), FUN = function(x){
        cod_pos = which(codons == x)
        if(length(cod_pos) == 0) {
            0
        } else {
            cod_pos
        }
    }))
    result
} 

simulate_subs = function(seqA, brlen, omega) {
    # encode sequence
    enc_seq = encode_seq(seqA)
    
    # deal with stop codons
    stops = which(enc_seq == 0)
    stopifnot(length(stops) <= 1)  # more than one stop codons
    if(length(stops) == 1) {
        stopifnot(stops == length(enc_seq))  # stop codon is not at the end
        stop_cod = substr(seqA, length(enc_seq) * 3 - 2, length(enc_seq) * 3)
        enc_seq = enc_seq[-length(enc_seq)]
    }

    # Evolve
    P = MG94(brlen, omega)
    
    codons = get_codonstr(61)
    evolved_seq = unlist(lapply(X = enc_seq, FUN = function(x) {
        sample(codons, size = 1, prob = P[x, ])
    }))
    
    # add ending stop codons
    if(length(stops) == 1) {
        evolved_seq = c(evolved_seq, stop_cod)
    }
    # convert back to string
    evolved_seq = str_c(evolved_seq, collapse = "")
    stopifnot(nchar(seqA) == nchar(evolved_seq))
    evolved_seq
}

simulate_indels = function(seqA, g, e) {
    # Triplet Model
    # g: gap open
    # e: gap extend
    
    # transition probabilities
    m2m = (1-g)^2
    m2i = g
    m2d = (1-g)*g
    i2m = (1-e)*(1-g)
    i2i = e
    i2d = (1-e)*g
    d2m = (1-e)
    d2i = 0
    d2d = e
    
    # transition probability matrix
    transitionP = matrix(data = c(m2m, m2i, m2d, i2m, i2i, i2d, d2m, d2i, d2d),
                         nrow = 3,
                         ncol = 3,
                         byrow = TRUE)
    
    # states
    states = c(M, I, D)
    names(states) = c("M", "I", "D")
    
    repeat{
        n_char = 0
        prev_state = 1
    
        path = c()
        while(n_char < nchar(seqA)) {
            current_state = sample(states, size = 1, prob = transitionP[prev_state, ])
            path[length(path)+1] = current_state
            prev_state = current_state
            if(current_state != I) {  # if current state is not insertion
                n_char = n_char + 1
            }
        }
        # if pattern doesn't have a match end transition redo
        if(sample(states, size = 1, prob = transitionP[prev_state, ]) != 1) {
            next
        }
        
        rle_res = rle(path)
        
        # if pattern doesn't start with a match redo
        if(rle_res$values[1] != M) {
            next
        }
        
        # if length of indels is not multiple of 3 redo
        skip = FALSE
        mapply(FUN = function(state, len) {
            if(state != M && len %%3 != 0) {
                skip <<- TRUE
            }
        }, rle_res$values, rle_res$lengths)
    
        if(!skip) {
            break
        }
    }
    
    rle_res
}
    
insert_cigar = function(rle_res, seqA, evolved_seq, omega) {
    nucs = c('A', 'C', 'G', 'T')
    p_nucs = c(0.308, 0.185, 0.199, 0.308)  # nucleotide stationary frequencies A,C,G,T
    
    A = seqA
    B = evolved_seq
    
    pos = cumsum(rle_res$lengths)-rle_res$lengths+1
    for(i in seq_along(rle_res$lengths)) {
        if(rle_res$values[i] == D) {  # deletion
            substr(B, pos[i], pos[i] + rle_res$lengths[i]) = paste0(rep('-', rle_res$lengths[i]), collapse = "")
        } else if(rle_res$values[i] == I) {  # insertion
            A = paste0(substr(A, 1, pos[i] - 1),
                       paste0(rep('-', rle_res$lengths[i]), collapse = ""),
                       substr(A, pos[i], nchar(A)))
            ins = sample(nucs, rle_res$lengths[i], prob = p_nucs, replace = TRUE)
            B = paste0(substr(B, 1, pos[i] - 1),
                       paste0(ins, collapse = ""),
                       substr(B, pos[i], nchar(B)))
        }
    }
    
    # check result
    stopifnot(nchar(A) == nchar(B))
    am = strsplit(A, split = "")[[1]] != "-"
    bm = strsplit(B, split = "")[[1]] != '-'
    h = ifelse(am, ifelse(bm, M, D), ifelse(bm, I, 0))
    h = rle(h)
    stopifnot(h$values == rle_res$values, h$lengths == rle_res$lengths)
    if(pass_selection(A, B, omega)) {
        return(as.list(c(A, B)))
    } else {
        return(list())
    }
}

check_selection = function(A, B, omega) {
    a = strsplit(A, split = "")[[1]]
    b = strsplit(B, split = "")[[1]]

    # list of synonymoyse codons
    codonstr = get_codonstr()
    syn = syncodons(codonstr)
    names(syn) = toupper(names(syn))
    syn = lapply(syn,toupper)

    # check gaps
    g = which(a == '-')
    if(length(g) == 0) {  # if not gaps return TRUE
        return(TRUE)
    }
    g.dif = diff(g)
    g.breaks = which(g.dif != 1)
    g.breaks = c(g.breaks, length(g))
    pos = g[1]
    for(i in g.breaks) {
        if(pos %% 3 == 1) {  # if phase 0
            if(runif(1) < omega) {  # check for type I (syn)
                return(FALSE)
            }
            pos = g[i + 1]
            next
        }

        # 123 123
        # grab complete codons around phase1 and phase2
        if(pos %% 3 == 2) {  # phase1
            p.start = pos - 1
            p.end = g[i] + 2
        } else {  # phase2
            p.start = pos - 2
            p.end = g[i] + 1
        }
        codA = a[p.start:p.end]
        codA = paste0(codA[codA != '-'], collapse = "")
        codB.char = b[p.start:p.end]
        codB = c()
        for(i in seq(1, length(codB.char), 3)) {
            codB = c(codB, paste0(codB.char[i:(i+2)], collapse = ""))
        }
        if(any(codB %in% syn[codA][[1]])) {  # check for type I (syn)
            if(runif(1) < omega) {
                return(FALSE)
            }
        } else {
            if(runif(1) > omega) {  # check for type II (nonsyn)
                return(FALSE)
            }
        }
        pos = g[i + 1]
    }
    return(TRUE)
}

pass_selection = function(A, B, omega) {
    if(!check_selection(A, B, omega)) {  # check insertions
        return(FALSE)
    }

    if(!check_selection(B, A, omega)) {  # check deletions
        return(FALSE)
    }

    return(TRUE)
}


simulate_triplet = function(input, output, ...) {
    # Triplet model substitution probabilities
    params = data.frame(brlen = 0.184, omega = 0.2, g = 0.001, e = 1 - 1/6)
    extra_args = list(...)
    if(length(extra_args) > 0) {
        for(i in length(extra_args)) {
            params[i] = as.numeric(extra_args[i])
        }
    }

    # Read sequence and process
    seqA = read.fasta(file = input,
                          as.string = TRUE,
                          seqonly = TRUE,
                          forceDNAtolower = FALSE)[[1]]
    
    set.seed(sum(utf8ToInt(basename(output))))
    evolved_seq = simulate_subs(seqA, params$brlen, params$omega)
    repeat{
        rle_res = simulate_indels(seqA, params$g, params$e)
        seqs = insert_cigar(rle_res, seqA, evolved_seq, params$omega)
        if(length(seqs) == 2) {
            break
        }
    }
    
    # fix issues with seqinr
    source("scripts/write_fasta.R")
    
    # save sequences
    write.fasta(sequences = seqs,
                seqnames = c("A", "B"),
                as.string = TRUE,
                file.out = output)
    
    # output to stdout CIGAR string and input file with original (ancestor) sequence
    cigar = paste0(rle_res$lengths, ifelse(rle_res$values == M, 'M', ifelse(rle_res$values == I, 'I', 'D')))
    print(paste0(
        basename(input), ",",
        paste0(cigar, collapse = ""), ",",
        params$brlen, ",",
        params$omega
    ))
}

if(!interactive()) {
	M = 1
	I = 2
	D = 3
    ARGS = commandArgs(trailingOnly = TRUE)
    len_args = length(ARGS)
    # args required: input, output
    # args optional: brlen, omega, g, e
    if(len_args == 2) {
        simulate_triplet(ARGS[1], ARGS[2])
    } else if(len_args > 2 && len_args <= 6) {
        simulate_triplet(ARGS[1], ARGS[2], ARGS[-c(1,2)])
    } else {
        stop("Invalid number of arguments.")
    }
}

