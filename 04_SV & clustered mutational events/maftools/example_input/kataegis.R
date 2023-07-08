#!/share/public/software/R/3.5.3/lib64/R/bin/Rscript
library(maftools)
args = commandArgs(trailingOnly=TRUE)

maf = read.maf(maf = args[1])
sample = args[2]


	out = rainfallPlot(maf = maf, tsb = sample, detectChangePoints = TRUE, pointSize = 0.4, savePlot=TRUE)
	save(out, file=paste0(sample,".rainfallPlot.Rdata"))

