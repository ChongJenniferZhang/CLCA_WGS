library("dndscv")
library(Biostrings)
args = commandArgs(trailingOnly=TRUE)
mutations = read.table(args[1] , header=T, colClasses=c("character","character","numeric","character","character"))
dndsout = dndscv(mutations)

write.table(dndsout$sel_cv, file="sel_cv.xls", sep="\t", col.names=T, row.names=F, quote=F)
#write.table(dndsout$globaldnds, file=globaldnds.xls, sep="\t", col.names=T, row.names=F, quote=F)
#write.table(dndsout$annotmuts, file=annotmuts.xls, sep="\t", col.names=T, row.names=F, quote=F)
#print(dndsout$nbreg$theta)

