#!/usr/bin/env Rscript

# set filename and codon
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("Usage: Rscript 03_plots.R statistics_file.txt [codon]")
}
fname <- args[1]
codon <- if (length(args) == 2) args[2] else "CTG"  # codon defaults to CTG

# read file
conn <- file(fname, open="rt")
x <- TRUE
while (x) {
    l <- readLines(conn, n=1)
    x <- !grepl("Translation", l)
}
df <- read.table(conn, sep="\t", header=FALSE)
colnames(df) <- unlist(strsplit(l, "\t"))
rownames(df) <- df[,1]
df[,1] <- NULL
close(conn)

# convert values to fraction
mat <- prop.table(data.matrix(df), 2)
cuts <- cut(mat[,1],
    breaks = seq(0, 1, by=0.1),
    include.lowest=TRUE)

# set color palette to gray, ranging from light to dark gray
cols <- gray.colors(10, start=1, end=0)
palette(cols)

# plot
outfile <- gsub(".txt", ".pdf", fname)
pdf(outfile)
barplot(mat[,1],
    col=as.factor(cuts),
    ylim=c(0,1),
    las=2,
    xlab="Translation",
    ylab="Fraction of PSMs",
    main=paste(codon, "translation"))
legend("topright",
    legend=c(1,NA,NA,NA,NA,0.5,NA,NA,NA,0),  # same length as cols
    fill=rev(cols),
    y.intersp=0.5,
    border=NA,
    bty="n")
dev.off()
