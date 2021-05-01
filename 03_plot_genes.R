#!/usr/bin/env Rscript

# set filename and codon
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript 03_plot_genes.R genes_file.csv [codon]")
}
fname <- args[1]
codon <- if (length(args) == 2) args[2] else "CTG"  # codon defaults to CTG

# read file
df <- read.csv(fname, check.names=FALSE)
sorteddf <- df[order(-df$`Sequence coverage [%]`),]

# plot
outfile <- gsub(".csv", ".pdf", fname)
pdf(outfile)
par(mfrow=c(5,1), mar=c(3,4,1,1), oma=c(2,2,0,0), xpd=NA)
barplot(sorteddf$`Sequence coverage [%]`,
        main="Sequence coverage",
        ylab="Percentage",
        ylim=c(0,100))
barplot(sorteddf$`#PSMs`,
        main="PSMs",
        ylab="Occurrence")
barplot(sorteddf$`#non-redundant peptides`,
        main="Non-redundant peptides",
        ylab="Occurrence")
barplot(sorteddf$`CTG coverage [%]`,
        main="CTG coverage",
        ylab="Percentage",
        ylim=c(0,100))
barplot(sorteddf$`#PSMs containing CTG`,
        main="PSMs containing CTG",
        ylab="Occurrence")
title(xlab="Gene index", line=1)
dev.off()
