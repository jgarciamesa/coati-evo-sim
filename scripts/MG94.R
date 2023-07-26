library(Matrix)
library(stringr)
library(seqinr)


get_codons = function(n = 64) {
    stopifnot(n == 61 || n == 64)
    nucs = c("A","C","G","T")
    codons = cbind(rep(nucs,each=16),
                   rep(nucs,times=4,each=4),
                   rep(nucs,16))
    if(n == 61) {
        codons = codons[-c(57, 51, 49), ]
    }
    return(codons)
}

get_codonstr = function(n = 64) {
    codons = get_codons(n)
    apply(codons, 1, str_c, collapse = "")
}

# original code from mutation.R - toycoati - Reed Cartwright
MG94 = function(brlen = 0.0133, omega = 0.2){
  # parameters
  # Yang (1994) Estimating the pattern of nucleotide substitution
  nucs = c("A","C","G","T")

  nuc_freqs = c(0.308,0.185,0.199,0.308)

  nuc_q = c(-0.818, 0.132, 0.586, 0.1,
            0.221, -1.349, 0.231, 0.897,
            0.909, 0.215, -1.322, 0.198,
            0.1, 0.537, 0.128, -0.765)
  nuc_q = matrix(nuc_q,4,4,byrow=T)

  #omega = 0.2
  #brlen = 0.0133


  # construct codons and figure out which ones are synonymous
  codons = get_codons()
  codonstrs = apply(codons,1,str_c,collapse="")
  codonstrs = codonstrs[-c(57, 51, 49)]

  syn = syncodons(codonstrs)
  names(syn) = toupper(names(syn))
  syn = lapply(syn,toupper)

  # MG94 model - doi:10.1534/genetics.108.092254
  Q = matrix(0,64,64)
  Pi = rep(0,64)
  # construct the transition matrix
  for(i in 1:64) {
    Pi[i] = prod(nuc_freqs[match(codons[i,],nucs)])
    for(j in 1:64) {
      if(i == j) {
        Q[i,j] = 0
      } else if(i %in% c(49, 51, 57) || j %in% c(49, 51, 57)) {
          next
      } else if(sum(codons[i,] != codons[j,]) > 1) {
        Q[i,j] = 0
      } else {
        if(codonstrs[j] %in% syn[[ codonstrs[i] ]]) {
          w = 1
        } else {
          w = omega
        }
        o = which(codons[i,] != codons[j,])
        x = which(nucs == codons[i,o])
        y = which(nucs == codons[j,o])

        Q[i,j] = w*nuc_q[x,y];
      }
    }
  }
  
  Q = Q[-c(57,51,49),]
  Q = Q[, -c(57,51,49)]
  Pi = Pi[-c(57,51,49)]

  # normalize Q
  diag(Q) = -rowSums(Q)
  Q = Q / -sum(Pi*diag(Q))

  # construct transition matrix
  P = expm(Q*brlen)
  return(P)
}

marMG94 = function(brlen = 0.0133, omega = 0.2, strategy = "sum") {
    P = MG94(brlen, omega)
    codons = get_codons()
    codons = codons[-c(57,51,49), ]
    nucs = c("A","C","G","T")

    p_marg = array(0, dim = c(61,3,4))
    if(strategy == "best" || strategy == "max") {
        for(i in 1:61) {
            for(j in 1:3) {
                for(k in 1:4) {
                    p_marg[i,j,k] = max(P[i,codons[,j]==nucs[k]])
                }
            }
        }
    } else { # sum
        for(i in 1:61) {
            for(j in 1:3) {
                for(k in 1:4) {
                    p_marg[i,j,k] = sum(P[i,codons[,j]==nucs[k]])
                }
            }
        }
    }
    p_marg
}
