#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args) < 2) {
    cat("USAGE: Nchart.R reads.lens out.pdf \n")
} else {

lensfile <- args[[1]]
outfile <- args[[2]]

cat("loading read lengths from: ", lensfile, "\n")
s <- scan(lensfile)

cat("creating plot in: ", outfile, "\n")

pdf(outfile, width=10, height=5)

a<-hist(s, breaks=500, xlim=c(0,50000), xlab="read length", ylab="Frequency", main="Read Length Histogram")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
abline(v=c(0,10000,20000,30000,40000,50000), col="white")
lines(a, col="blue")

dev.off()
}
