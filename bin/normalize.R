#!/usr/bin/Rscript
library("optparse")
library("mgcv")

option_list <- list(
         make_option(c("-i", "--infile"), action = "store", type = "character", help = "sambamba depth file name"),
         make_option(c("-p", "--prefix"), action = "store", type = "character", help = "output file prefix"),
         make_option(c("-d", "--datadir"), action = "store", type = "character", help = "data file dir")
         )
opt <- parse_args(OptionParser(option_list = option_list))

over_all_bed <- read.delim(paste0(opt$datadir, "/ct_panel_152_auto_withg_uniq.data"), header=F,col.names=c("X..chrom","chromStart","chromEnd","overlap_score","gc","round_gc","mapp"))

probe_stat <- read.delim(paste0(opt$datadir, "/probe_level_stat"), header=F, col.names=c("X..chrom","chromStart","chromEnd","gene","mean","stdev"))
probe_stat <- probe_stat[,-4]

plot_data <- read.delim(paste0(opt$datadir, "/plot.data"), header=F, col.names=c("gene","boundary_a","boundary_b","x","y"))

chr_data <-  read.delim(paste0(opt$datadir, "/chr.data"), header=F, col.names=c("chr","boundary","x"))


data <- read.delim(opt$infile, header = TRUE)
data$log_depth <- log(1 + data$meanCoverage)

newdata <-merge(data,over_all_bed,by=c("X..chrom","chromStart","chromEnd"))

model <- gam(log_depth ~ s(gc,bs='ps') + s(round_gc,bs='ps') + s(overlap_score,bs='tp') + mapp , data=newdata)

predict_depth <- exp(predict(model,newdata)) - 1

newdata$logratio <- log2(newdata$meanCoverage / predict_depth)

###use a qunatile(0.2,0.3) to compensate the mean
offset = quantile(newdata$logratio, probs=seq(0.2,0.3,0.01))
newdata$logratio = newdata$logratio - mean(offset)

newdata <- newdata[order(as.numeric(substr(newdata$X..chrom,4,10)),newdata$chromStart),]

plotdata <- merge(newdata,probe_stat,by=c("X..chrom","chromStart","chromEnd"),all.x=T)
plotdata[is.na(plotdata)] <- 0
plotdata <- plotdata[order(as.numeric(substr(plotdata$X..chrom,4,10)),plotdata$chromStart),]

prefix = opt$prefix
write.table(newdata,paste(prefix,".tsv",sep=""),sep="\t",quote=F,row.names=F)
write.table(plotdata,paste(prefix,".plot.tsv",sep=""),sep="\t",quote=F,row.names=F)
png(paste(prefix,".png",sep=""),1500,900)
plot(seq(length(plotdata$logratio)),plotdata$logratio - plotdata$mean,xlab="probe position",ylab="normalized cov",ylim=c(-1.35,3.6))
abline(v=c(plot_data$boundary_a,plot_data$boundary_b), col="green")
text(x=plot_data$x, y=plot_data$y, labels=plot_data$gene, pos=1, adj=0, col="red")
abline(v=chr_data$boundary, col="blue")
text(x=chr_data$x, labels=chr_data$chr, y=3.6, pos=1, adj=0, col="blue")
dev.off()
