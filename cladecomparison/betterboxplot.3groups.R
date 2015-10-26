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
	# if there is an outlier, label it
	for (i in 1:length(set1)) {
		if(!is.na(set1[i])) {
			if (abs(as.numeric(set1[i])-median(set1, na.rm=T))>2*sd(set1, na.rm=T)) {
				text(x=0.9, y=set1[i], pos=1, labels=strsplit(names(set1)[i], "_")[[1]][1])
			}
		}
	}
	for (i in 1:length(set2)) {
		if(!is.na(set2[i])) {
			if (abs(as.numeric(set2[i])-median(set2, na.rm=T))>2*sd(set2, na.rm=T)) {
				text(x=1.9, y=set2[i], pos=1, labels=strsplit(names(set2)[i], "_")[[1]][1])
			}
		}
	}
	for (i in 1:length(set3)) {
		if(!is.na(set3[i])) {
			if (abs(as.numeric(set3[i])-median(set3, na.rm=T))>2*sd(set3, na.rm=T)) {
				text(x=2.9, y=set3[i], pos=1, labels=strsplit(names(set3)[i], "_")[[1]][1])
			}
		}
	}
	dev.off()
}

enviro <- read.delim("path-env.NAs.genes", header=F, stringsAsFactors=F)
rhizo <- read.delim("path-rhiz.NAs.genes", header=F, stringsAsFactors=F)
genes <- rbind(enviro, rhizo)
genes <- genes[!duplicated(genes[,1]),]

#for (i in 1:nrow(genes)) {
#	gene <- genes[i,1]
#	name <- genes[i,2]
#	print(gene)
#	plotbox(gene, name)
#}


plotposterbox <- function(gene, name) {
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
	pdf(paste("boxplots/", type, "/", gene, ".pdf", sep=""), width=5, height=5)
	plot(x=seq(0.9,1.1,0.2/(length(set1)-1)), y=set1, xlim=c(0.7,2.3), ylim=c(min(pooled, na.rm=T)-25, max(pooled, na.rm=T)+25), main=paste("Score distributions for\n", name, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.5, cex.lab=1.3, cex.main=1.5, cex.axis=1.3, frame=F, col="coral4")
	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, col="darkolivegreen", pch=16, cex=1.5)
	points(x=seq(1.9,2.1,0.2/(length(set3)-1)), y=set3, col="darkgoldenrod4", pch=16, cex=1.5)
	#	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, col="darkolivegreen", pch=16, cex=1.3)
	#points(x=seq(2.9,3.1,0.2/(length(set3)-1)), y=set3, col="darkolivegreen", pch=16, cex=1.3)
	axis(side=1, at=c(1,2), labels=c("Pathogenic", "Non-pathogenic"), cex.axis=1.3)
	polygon(x=c(0.8, 0.8, 1.2, 1.2), y=c(rep(median(set1, na.rm=T), 4)), col="black", lwd=3)
	polygon(x=c(1.8, 1.8, 2.2, 2.2), y=c(rep(median(c(set2, set3), na.rm=T), 4)), col="black", lwd=3)
	#polygon(x=c(2.8, 2.8, 3.2, 3.2), y=c(rep(median(set3, na.rm=T), 4)), col="black", lwd=4)
	dev.off()
	png(paste("boxplots/", type, "/", gene, ".png", sep=""), width=400, height=400)
	plot(x=seq(0.9,1.1,0.2/(length(set1)-1)), y=set1, xlim=c(0.7,2.3), ylim=c(min(pooled, na.rm=T)-25, max(pooled, na.rm=T)+25), main=paste("Score distributions for\n", name, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.5, cex.lab=1.3, cex.main=1.5, cex.axis=1.3, frame=F, col="coral4")
	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, col="darkolivegreen", pch=16, cex=1.5)
	points(x=seq(1.9,2.1,0.2/(length(set3)-1)), y=set3, col="darkgoldenrod4", pch=16, cex=1.5)
	#	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, col="darkolivegreen", pch=16, cex=1.3)
	#points(x=seq(2.9,3.1,0.2/(length(set3)-1)), y=set3, col="darkolivegreen", pch=16, cex=1.3)
	axis(side=1, at=c(1,2), labels=c("Pathogenic", "Non-pathogenic"), cex.axis=1.3)
	polygon(x=c(0.8, 0.8, 1.2, 1.2), y=c(rep(median(set1, na.rm=T), 4)), col="black", lwd=3)
	polygon(x=c(1.8, 1.8, 2.2, 2.2), y=c(rep(median(c(set2, set3), na.rm=T), 4)), col="black", lwd=3)
	#polygon(x=c(2.8, 2.8, 3.2, 3.2), y=c(rep(median(set3, na.rm=T), 4)), col="black", lwd=4)
	dev.off()
}

for (i in 1:nrow(genes)) {
	gene <- genes[i,1]
	name <- genes[i,2]
	print(gene)
	plotposterbox(gene, name)
}

# path-nonpath analysis

genes <- read.delim("path-nonpath.noDBS.genes", header=F, stringsAsFactors=F)

plotposterbox <- function(gene, name) {
	set1 <- pathscores2[,gene]
	set2 <- rep(NA, nrow(rhizscores2))
	set3 <- rep(NA, nrow(envscores2))
	if(!is.na(match(gene, pathrhiz[,1]))) {
		if (!is.na(match(pathrhiz[pathrhiz[,1]==gene,2], colnames(rhizscores2)))) {
			set2 <- rhizscores2[,pathrhiz[pathrhiz[,1]==gene,2]]
			if (!is.na(match(gene, rhizo[,1]))) {
			}
		}
	}
	if (!is.na(match(gene, pathenv[,1]))) {
		if(!is.na(match(pathenv[pathenv[,1]==gene,2], colnames(envscores2)))) {
			set3 <- envscores2[,pathenv[pathenv[,1]==gene,2]]
			if (!is.na(match(gene, enviro[,1]))) {
			}
		}
	}
#	set1 <- g1s[,gene]
#	set2 <- g2s[,pathenv[pathenv[,1]==gene,2]]
	pooled <- append(set1, set2)
	pooled <- append(pooled, set3)
	pdf(paste("boxplots/noDBS/", gene, ".pdf", sep=""), width=5, height=5)
	plot(x=seq(0.9,1.1,0.2/(length(set1)-1)), y=set1, xlim=c(0.7,2.3), ylim=c(min(pooled, na.rm=T), max(pooled, na.rm=T)), main=paste("Score distributions for\n", name, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.5, cex.lab=1.3, cex.main=1.2, cex.axis=1.3, frame=F, col="coral4")
	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, col="darkolivegreen", pch=16, cex=1.5)
	points(x=seq(1.9,2.1,0.2/(length(set3)-1)), y=set3, col="darkgoldenrod4", pch=16, cex=1.5)
	axis(side=1, at=c(1,2), labels=c("Pathogenic", "Non-pathogenic"), cex.axis=1.3)
	polygon(x=c(0.8, 0.8, 1.2, 1.2), y=c(rep(median(set1, na.rm=T), 4)), col="black", lwd=3)
	polygon(x=c(1.8, 1.8, 2.2, 2.2), y=c(rep(median(c(set2, set3), na.rm=T), 4)), col="black", lwd=3)
	dev.off()
	png(paste("boxplots/noDBS/", gene, ".png", sep=""), width=400, height=400)
	plot(x=seq(0.9,1.1,0.2/(length(set1)-1)), y=set1, xlim=c(0.7,2.3), ylim=c(min(pooled, na.rm=T), max(pooled, na.rm=T)), main=paste("Score distributions for\n", name, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.5, cex.lab=1.3, cex.main=1.2, cex.axis=1.3, frame=F, col="coral4")
	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, col="darkolivegreen", pch=16, cex=1.5)
	points(x=seq(1.9,2.1,0.2/(length(set3)-1)), y=set3, col="darkgoldenrod4", pch=16, cex=1.5)
	axis(side=1, at=c(1,2), labels=c("Pathogenic", "Non-pathogenic"), cex.axis=1.3)
	polygon(x=c(0.8, 0.8, 1.2, 1.2), y=c(rep(median(set1, na.rm=T), 4)), col="black", lwd=3)
	polygon(x=c(1.8, 1.8, 2.2, 2.2), y=c(rep(median(c(set2, set3), na.rm=T), 4)), col="black", lwd=3)
	dev.off()
}

for (i in 1:nrow(genes)) {
	gene <- genes[i,1]
	name <- genes[i,2]
	print(gene)
	plotposterbox(gene, name)
}

