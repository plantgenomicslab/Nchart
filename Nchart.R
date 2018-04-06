#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args) < 2) {
  cat("USAGE: Nchart.R out.png lens1 lens2 lens3 ... lensN\n")
} else {

# output prefix
outfile <- args[[1]]

query_filenames = tail(args, -1)
query_names <- query_filenames

## Assume the first entry is the reference
filename_ref <- query_filenames[[1]]

cat("Creating Nchart\n")
cat("===============\n")
cat("outfile: ",         outfile, "\n")
cat("filename_ref: ",    filename_ref, "\n")
cat("query_filenames: ", query_filenames, "\n")
cat("query_names: ",     query_names, "\n")


stop()

library(ggplot2)
library(scales)


# first color is the reference, the others are the queries in order
colors <- rainbow(length(query_filenames))



#############################################################################

ref.data <- read.table(filename_ref, header=FALSE)
names(ref.data) <- c("name","length")
ref.data$length <- sort(as.numeric(ref.data$length),decreasing=TRUE)

genome.length <- max(sum(ref.data$length))


ref.cumsum <- data.frame(NG=cumsum(ref.data$length/genome.length*100),contig.length=ref.data$length,contig.source="Reference")

ref.cumsum.0 <- rbind(data.frame(NG=c(0),contig.length=max(ref.cumsum$contig.length),contig.source="Reference"),ref.cumsum)


for (i in seq(length(query_filenames))) {
  filename_query = query_filenames[i]
  name_query = query_names[i]
  print(paste(filename_query,name_query))
  
  query.data <- read.table(filename_query, header=FALSE)
  
  names(query.data) <- c("name","length")
  
  query.data$length <- sort(as.numeric(query.data$length),decreasing=TRUE)
  
  query.cumsum <- data.frame(NG=cumsum(query.data$length/genome.length*100),contig.length=query.data$length,contig.source=name_query)
  
  query.cumsum.0 <- rbind(data.frame(NG=c(0),contig.length=max(query.cumsum$contig.length),contig.source=name_query),query.cumsum)
  
  ref.cumsum <- rbind(ref.cumsum,query.cumsum)
  ref.cumsum.0 <- rbind(ref.cumsum.0,query.cumsum.0)
}

bp_format<-function(num) {
  if (num > 1000000000) {
    paste(formatC(num/1000000000,format="f",digits=1,big.mark=",",drop0trailing = TRUE)," Gbp",sep="")
  }
  else if (num > 1000000) {
    paste(formatC(num/1000000,format="f",digits=1,big.mark=",",drop0trailing = TRUE)," Mbp",sep="")
  }
  else {
    paste(formatC(num,format="f",big.mark=",",drop0trailing = TRUE), " bp", sep="")
  } 
}

theme_set(theme_bw(base_size = 12) + theme(panel.grid.minor = element_line(colour = NA)))


percent_format <- function(num) {
  paste(num,"%", sep="")
}

png(file=outfile, width=1600,height=1000,res=200)

print(
  ggplot(ref.cumsum.0, aes(x = NG, y = contig.length, color=contig.source)) + 
    #             xlim(0,100) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits=c(1,genome.length)) + 
    scale_x_continuous(labels=percent_format) +
    geom_path(size=1.5,alpha=0.5) + 
    geom_point(data=ref.cumsum, size=2,alpha=0.5) + 
    labs(x = paste("Percentage of reference (",bp_format(genome.length),")", sep=""),y="Sequence length",colour="Assembly",title="Cumulative sequence length") +
    scale_color_manual(values=colors) +
    annotation_logticks(sides="lr")
)

dev.off()

}


