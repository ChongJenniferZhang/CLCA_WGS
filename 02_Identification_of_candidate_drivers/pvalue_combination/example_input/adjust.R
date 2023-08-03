args=commandArgs(TRUE)
d<-read.table(args[1], header=T)
d=d[order(d$Brown_observed),]
p=d$Brown_observed
adj = p.adjust(p, method='BH')
d$q = adj

write.table(d,file=args[2], quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
