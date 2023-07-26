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
