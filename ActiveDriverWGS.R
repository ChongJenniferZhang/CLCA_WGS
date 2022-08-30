#!/share/public/software/R/3.5.3/lib64/R/bin/Rscript
library("ActiveDriverWGS")
args = commandArgs(trailingOnly=TRUE)

mutations = read.table("mutations" , header=T, colClasses=c("character","numeric","numeric","character","character","character"))
elements = read.table("elements" , header=T, colClasses=c("character","numeric","numeric","character"))
result = ActiveDriverWGS(mutations, elements, sites = NULL, window_size = 50000,
					  		filter_hyper_MB = 30, recovery.dir = "recovery_one", mc.cores = 10)

write.table(result, file=args[1], sep="\t", col.names=T, row.names=F, quote=F)
