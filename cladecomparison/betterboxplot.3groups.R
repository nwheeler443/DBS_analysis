# give a gene name and a boxplot will be generated that can be used for presentation purposes
# some genes are on orthlist but not dbs because of incompatible domain architectures

pathenv <- read.delim("path-env.dbs/orthlist.dbs", header=F, stringsAsFactors=F)
pathrhiz <- read.delim("path-rhiz.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("Pathogenic.scores.NAs.txt", header=T, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
envscores <- read.delim("Environmental.scores.NAs.txt", header=T, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)
rhizscores <- read.delim("Rhizosphere.scores.NAs.txt", header=T, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

set1 <- vector()
set2 <- vector()
set3 <- vector()

plotbox <- function(gene, name) {
	set1 <- pathscores2[,gene]
	set2 <- rep(NA, nrow(rhizscores2))
	set3 <- rep(NA, nrow(envscores2))
	type = ""
	if(!is.na(match(gene, pathrhiz[,1]))) {
		if (!is.na(match(pathrhiz[pathrhiz[,1]==gene,2], colnames(rhizscores2)))) {
			set2 <- rhizscores2[,pathrhiz[pathrhiz[,1]==gene,2]]
			if (!is.na(match(gene, rhizo[,1]))) {
				type = "path-rhiz"
			}
		}
	}
	if (!is.na(match(gene, pathenv[,1]))) {
		if(!is.na(match(pathenv[pathenv[,1]==gene,2], colnames(envscores2)))) {
			set3 <- envscores2[,pathenv[pathenv[,1]==gene,2]]
			if (!is.na(match(gene, enviro[,1]))) {
				type = "path-env"
			}
		}
	}
	if (!is.na(match(gene, rhizo[,1])) & !is.na(match(gene, enviro[,1]))) {
		type = "shared"
	}
	pooled <- append(set1, set2)
	pooled <- append(pooled, set3)
	png(paste("boxplots/", type, "/", gene, ".png", sep=""))
	plot(x=c(seq(0.9,1.1,0.2/(length(set1)-1)), seq(1.9,2.1,0.2/(length(set2)-1)), seq(2.9,3.1,0.2/(length(set3)-1))), y=c(set1, set2, set3), xlim=c(0.5,3.5), ylim=c(min(pooled, na.rm=T)-25, max(pooled, na.rm=T)+25), main=paste("Score distributions for gene: \n", name, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.3, cex.lab=1.3, cex.main=1.5, cex.axis=1.3, frame=F)
	axis(side=1, at=c(1,2,3), labels=c("Pathogenic", "Rhizosphere", "Environmental"), cex.axis=1.3)
	polygon(x=c(0.8, 0.8, 1.2, 1.2), y=c(rep(median(set1, na.rm=T), 4)), col="black", lwd=4)
	polygon(x=c(1.8, 1.8, 2.2, 2.2), y=c(rep(median(set2, na.rm=T), 4)), col="black", lwd=4)
	polygon(x=c(2.8, 2.8, 3.2, 3.2), y=c(rep(median(set3, na.rm=T), 4)), col="black", lwd=4)
	dev.off()
}

enviro <- read.delim("path-env.NAs.genes", header=F, stringsAsFactors=F)
rhizo <- read.delim("path-rhiz.NAs.genes", header=F, stringsAsFactors=F)
genes <- rbind(enviro, rhizo)
genes <- genes[!duplicated(genes[,1]),]
#genes <- first[!is.na(match(first[,1], others[,1])),]

for (i in 1:nrow(genes)) {
	gene <- genes[i,1]
	name <- genes[i,2]
	print(gene)
	plotbox(gene, name)
}
