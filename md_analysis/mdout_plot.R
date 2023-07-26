#!/blue/meilerlab/apps/Linux2/x86_64/R/2.15.0/bin/Rscript --vanilla
library(ggplot2,quietly=TRUE)

args <- commandArgs(TRUE)
name <- args[1]
data <-read.table(args[1], col.names=c("time", "prop"), fill=TRUE)

png(paste(args[1],"png",sep = "."))

ggplot(data, aes(time, prop)) + 
	geom_point() + 
    ggtitle(name) +
	theme_bw()

dev.off()
