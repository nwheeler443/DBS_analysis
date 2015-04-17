# give a gene name and a boxplot will be generated that can be used for presentation purposes

pathenv <- read.delim("path-env.dbs/orthlist.dbs", header=F, stringsAsFactors=F)
pathrhiz <- read.delim("path-rhiz.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("pathogenscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
envscores <- read.delim("environmentalscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)
rhizscores <- read.delim("rhizospherescores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

plotbox <- function(gene, name) {
	set1 <- pathscores2[,gene]
	if (length(match(pathenv[pathenv$V1==gene,2], colnames(envscores2)))>0) {
		if(length(match(pathrhiz[pathrhiz$V1==gene,2], colnames(rhizscores2)))>0) {
			set2 <- envscores2[,pathenv[pathenv[,1]==gene,2]]
			set3 <- rhizscores2[,pathrhiz[pathrhiz[,1]==gene,2]]
			pooled <- append(set1, set2, set3)
			png(paste("boxplots/Figures/shared/", gene, ".png", sep=""))
			plot(x=c(seq(0.9,1.1,0.2/(length(set1)-1)), seq(1.9,2.1,0.2/(length(set2)-1)), seq(2.9,3.1,0.2/(length(set3)-1))), y=c(set1, set2, set3), xlim=c(0.5,3.5), ylim=c(min(pooled, na.rm=T)-25, max(pooled, na.rm=T)+25), main=paste("Score distributions for gene \n", gene, "\n", name, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.3, cex.lab=1.3, cex.main=1.3, cex.axis=1.3, frame=F)
			axis(side=1, at=c(1,2,3), labels=c("Pathogenic", "Environmental", "Rhizosphere"), cex.axis=1.3)
			dev.off()
		}
		else {
			print("couldn't find gene 1\n")
		}
	}
	else {
		print("couldn't find genes\n")
	}
}

genes <- read.delim("path-env.NAs.genes", header=F, stringsAsFactors=F)
for (i in 1:nrow(genes)) {
	gene <- genes[i,1]
	name <- genes[i,2]
	plotbox(gene, name)
}
genes <- read.delim("path-rhiz.NAs.genes", header=F, stringsAsFactors=F)
for (i in 1:nrow(genes)) {
	gene <- genes[i,1]
	name <- genes[i,2]
	plotbox(gene, name)
}
