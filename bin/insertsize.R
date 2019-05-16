#!/usr/bin/Rscript

library("ggplot2")
library("optparse")

option_list <- list(
	make_option(c("-i", "--infile"), action = "store", type = "character", help = "samtools stat insert size file"),
	make_option(c("-p", "--prefix"), action = "store", type = "character", help = "output file prefix")
	)

opt <- parse_args(OptionParser(option_list = option_list))
data <- read.delim(opt$infile, header=F)

a <- sum(data[data$V1>=141&data$V1<=152,]$V3)
b <- sum(data[data$V1>=165&data$V1<=250,]$V3)
tumc <- ifelse(a/b < 0.21, 0.0123, 1.23*(a/b-0.2))
tumc <- ifelse(tumc > 0.9, 0.9, tumc)

plot <- ggplot(data=data, aes(x=V1,y=V3)) + geom_path() + geom_vline(data=data.frame(cv=c(122,134,141,152,165)),aes(xintercept=cv))
plot <- plot + geom_text(x=210,y=max(data$V3)*0.95,label= paste("tumc =",round(tumc,5)*100, "%"),size=9) + xlim(100,250)
plot <- plot + theme(legend.position="none")

png(paste(opt$prefix,"-dist.png",sep=""),800,360)
print(plot)
dev.off()

