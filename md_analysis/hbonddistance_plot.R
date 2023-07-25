#!/blue/meilerlab/apps/Linux2/x86_64/R/2.15.0/bin/Rscript --vanilla
library(Hmisc,quietly=TRUE)

args <- commandArgs(TRUE)
data <-read.table(args[1], fill=TRUE, header=TRUE)
head(data)

png(paste(args[3], args[2], "hd", "png",sep = "."))

hist.data.frame(data, main=colnames(data))
text(x = 1, y = 0, labels = args[2])

dev.off()