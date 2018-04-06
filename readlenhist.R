#s <- scan("data.txt")

pdf("skbr3.ont.pdf")

#pdf("skbr3.pacbio.pdf")
#a<-hist(s, breaks=500, xlim=c(0,50000))
a<-hist(s, breaks=2000, xlim=c(0,50000), xlab="read length", ylab="Frequency", main="Read Length Histogram")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
abline(v=c(0,10000,20000,30000,40000,50000), col="white")
lines(a, col="blue")

dev.off()
