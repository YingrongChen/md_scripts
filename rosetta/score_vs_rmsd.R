library(ggplot2,quietly=TRUE)

args <- commandArgs(TRUE)

data <-read.table(args[1],header=TRUE)

head(data)

png(paste(args[1],".png",sep = "_"))

ggplot(data, aes(total_score, rmsd)) + 
	geom_point() + 
	xlab("total_score") +
	theme_bw()

dev.off()
