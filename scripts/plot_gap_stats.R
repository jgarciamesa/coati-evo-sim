library(ggplot2)

plot_freq = function(brlen) {
    models = c("reference", "tri-mg", "mar-mg-sum", "mar-mg-max")
    freq = c()
    for(m in models) {
        f = read.csv(file = paste0("results/", brlen, "/gap_stats/freq-", m, ".csv"), header = TRUE)
        f = cbind(f, rep(m, nrow(f)))
        names(f) = c("gap_length", "count", "model")
        freq = rbind(freq, f)
    }
    
    pfreq = ggplot() +
        geom_point(data = freq, mapping = aes(x = gap_length, y = count, color = model))
    
    ggsave(filename = paste0("results/", brlen, "/figures/gap-freq.pdf"), plot = pfreq)
}

plot_phase = function(brlen) {
    models = c("reference", "tri-mg", "mar-mg-sum", "mar-mg-max")
    phase = c()
    for(m in models) {
        p = read.csv(file = paste0("results/", brlen, "/gap_stats/phase-", m, ".csv"), header = TRUE)
        p = cbind(p, rep(m, nrow(p)))
        names(p) = c("phase", "count", "model")
        phase = rbind(phase, p)
    }

    pphase = ggplot(data = phase, mapping = aes(x = model, y = count, fill = phase)) +
        geom_bar(stat = "identity", position = "stack")
    
    ggsave(filename = paste0("results/", brlen, "/figures/gap-phase.pdf"), plot = pphase)
}

plot_pos = function(brlen) {
    models = c("reference", "tri-mg", "mar-mg-sum", "mar-mg-max")
    pos = c()
    for(m in models) {
        p = read.csv(file = paste0("results/", brlen, "/gap_stats/pos-", m, ".csv"), header = TRUE)
        p = cbind(p, rep(m, nrow(p)))
        names(p) = c("position", "count", "model")
        pos = rbind(pos, p)
    }

    ppos = ggplot(data = pos, mapping = aes(x = position, y = count, color = model)) +
        geom_line(stat = "identity")
  
    ggsave(filename = paste0("results/", brlen, "/figures/gap-pos.pdf"), plot = ppos)
}

if(!interactive()) {
    ARGS = commandArgs(trailingOnly = TRUE)
    stopifnot(length(ARGS) == 2)
    if(ARGS[1] == "freq") {
        plot_freq(ARGS[2])
    } else if(ARGS[1] == "phase") {
        plot_phase(ARGS[2])
    } else if(ARGS[1] == "pos") {
        plot_pos(ARGS[2])
    } else {
        warning(paste("Option", ARGS[1], "not supported."))
    }
}