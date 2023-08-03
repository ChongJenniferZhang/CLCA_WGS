#!/share/public/software/R/3.5.3/lib64/R/bin/Rscript
library(ShatterSeek)
args = commandArgs(trailingOnly=TRUE)

SV <- read.table(args[1], header=T, sep="\t")
CN <- read.table(args[2], header=T, sep="\t")
sample=args[3]

SV_data <- SVs(chrom1=as.character(SV$chrom1),
	       pos1=as.numeric(SV$pos1),
	       chrom2=as.character(SV$chrom2),
	       pos2=as.numeric(SV$pos2),
	       SVtype=as.character(SV$SVtype),
	       strand1=as.character(SV$strand1),
	       strand2=as.character(SV$strand2))

CN_data <- CNVsegs(chrom=as.character(CN$chrom),
		   start=CN$start,
		   end=CN$end,
		   total_cn=CN$CN)

start_time <- Sys.time()
chromothripsis <- shatterseek(SV.sample=SV_data, seg.sample=CN_data)
save(chromothripsis, file=paste(sample,"chromothripsis.Rdata",sep="."))
end_time <- Sys.time()
print(paste0("Running time (s): ",round(end_time - start_time,digits=2)))

write.table(chromothripsis@chromSummary, file=paste(sample, "chromSummary.txt", sep="."), sep="\t", col.names=T, row.names=F, quote=F)

chromothripsis@chromSummary[!is.na(chromothripsis@chromSummary$start), ]

library(gridExtra)
library(cowplot)
library(ggplot2)

names=read.table("samplename.txt")
sample_name=as.character(names[names[,2]==sample,1])

for (ch in as.character(chromothripsis@chromSummary[!is.na(chromothripsis@chromSummary$start), ]$chrom)){
	plots_chr = plot_chromothripsis(ShatterSeek_output = chromothripsis,chr = ch, sample_name=sample)
	plot_chr = arrangeGrob(plots_chr[[1]], plots_chr[[2]], plots_chr[[3]], plots_chr[[4]], nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))
	ggsave(plot_chr, file=paste(sample,ch,"pdf", sep="."))
}




