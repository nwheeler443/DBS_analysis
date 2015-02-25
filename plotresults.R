plotresults <- function(comm, path, xaxvalues=c(-10,10), breaks=500) {
    dbs <- read.table(paste(comm, "-", path, "/results.dbs", sep=""))
    info <- read.delim(paste("genomes/", path, ".info", sep=""), header=F)
    
	genelist <- genelist <- cbind(dbs[,c(1,2,3,10,12)], info[match(dbs[,2], info[,1]), 2:ncol(info)])
    write.table(genelist, file=paste(comm, "-", path, "/genelist.txt", sep=""), quote=F, sep="\t", col.names=F, row.names=F)
    
    png(paste(comm, "-", path, "/graph.png", sep=""), width=1000, height=1000, res=150)
    hist(dbs[dbs$V10!=0,10], breaks=breaks, xlim=xaxvalues, main=paste("Results for ", comm, " vs ", path, sep=""), xlab="Delta bitscore")
    library(moments)
    text((xaxvalues[2]/2), nrow(dbs[-0.5<dbs$V10 & dbs$V10<0.5,])/2, paste("skewness = ", round(skewness(dbs[,10]), 2)))
    dev.off()
    
}