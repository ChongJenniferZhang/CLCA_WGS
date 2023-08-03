rm(list=ls()) 
args<-commandArgs(trailingOnly = T)

library("sequenza")
data.file<-args[1]
id<-args[2]
seqz.data <- read.seqz(data.file) 
str(seqz.data, vec.len = 2)

gc.stats <- gc.sample.stats(data.file)
str(gc.stats)
gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio /gc.vect[as.character(seqz.data$GC.percent)]
#pdf("GC_ration.pdf")
#par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
#matplot(gc.stats$gc.values, gc.stats$raw, type = 'b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2), xlab = 'GC content (%)', ylab = 'Uncorrected depth ratio')
#legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
#hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio,breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2),xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
#dev.off()

test <- sequenza.extract(data.file,assembly = "hg19")
names(test)
#chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]], ratio.windows = test$ratio[[1]], min.N.ratio = 1, segments = test$segments[[1]], main = test$chromosomes[1])
CP.example <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id=id)
