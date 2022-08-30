args = commandArgs(trailingOnly=TRUE)
library(MutationTimeR)
vcf <- readVcf(args[1], genome="GRCh37")
bb <- MutationTimeR:::loadBB(args[2])
purityPloidy <- read.table(args[3], header=TRUE, sep="\t")
sample <- args[4]
print(sample)
purityPloidyss = purityPloidy[purityPloidy$sample==sample, ]
purity <- purityPloidyss$purity
bb$clonal_frequency <- purity

gender='female'
if (purityPloidyss$gender == "MALE"){
	gender='male'
}
isWgd=FALSE
if (purityPloidyss$wholeGenomeDuplication == "true"){
	isWgd=TRUE
}

print(c(sample, gender, purity, isWgd))

mt <- mutationTime(vcf, bb, n.boot=1000, gender=gender, isWgd=isWgd)
save(mt, file=paste(sample, "mt.Rdata", sep="."))
print(1)
#write.table(mt$V, file=paste(sample,"V.txt",sep="."), sep="\t", col.names=T, row.names=F, quote=F)
# Classify as basic clonal states
write.table(table(mt$V$CLS), file=paste(sample,"classifyMutations.txt",sep="."), sep="\t", col.names=T, row.names=F, quote=F)

info(header(vcf)) <- rbind(info(header(vcf)),mtHeader())
info(vcf) <- cbind(info(vcf), mt$V)
MutationTimeR:::writeVcf(vcf, paste(sample, "info.vcf",sep="."))
print(2)

mcols(bb) <- cbind(mcols(bb),mt$T)
write.table(bb[,-5], file=paste(sample, "bb.txt",sep="."), sep="\t", col.names=T, row.names=F, quote=F)
print(3)

pdf(paste(sample, ".plot.pdf", sep="."))
plotSample(vcf,bb)
dev.off()
save(vcf, file=paste(sample, "vcf.Rdata", sep="."))
save(bb, file=paste(sample, "bb.Rdata", sep="."))
