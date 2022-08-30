library(sva)

d=read.table("raw.txt", header=T, row.names=1)

cdata <- as.matrix(d)

csif <- read.table("group.txt", header = T, sep = "\t", colClasses=c("character","numeric","numeric","numeric","character"))

modcombat = model.matrix(~cc, data = csif)

batch = csif$batch

#pdf("rmbatch.before.pdf", 20, 5)
#dist_mat <- dist(t(cdata))
#clustering <- hclust(dist_mat, method = "complete")
#plot(clustering, labels = batch, ylab="", xlab="", sub="")
#plot(clustering, labels = csif$label, ylab="", xlab="", sub="")
#dev.off()
#
#pdf("combat.pdf")
combat_edata = ComBat(dat=cdata+0.001, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
#write.table(combat_edata, file="combat_edata.xls", sep="\t", col.names=T, row.names=T, quote=F)
#dev.off()

#pdf("rmbatch.after.pdf", 20, 5)
#dist_mat_combat <- dist(t(combat_edata))
#clustering_combat <- hclust(dist_mat_combat, method = "complete")
#plot(clustering_combat, labels = batch, ylab="", xlab="", sub="")
#plot(clustering_combat, labels = csif$label, ylab="", xlab="", sub="")
#dev.off()

library("ggpubr")
library("ggthemes")
library("Rtsne")
d=combat_edata
set.seed(1)
tsne.info = Rtsne(t(d), perplexity=3)
colnames(tsne.info$Y) = c("tSNE_1","tSNE_2")
tSNE.data = data.frame(sample=colnames(d), Type=as.character(batch), tsne.info$Y)


pdf("tSNE_rmbatch.pdf")
ggscatter(tSNE.data, x="tSNE_1", y="tSNE_2", color="Type")
dev.off()
