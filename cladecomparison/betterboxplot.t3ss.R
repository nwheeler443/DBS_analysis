# give a gene name and a boxplot will be generated that can be used for presentation purposes

pathscores <- read.delim("T3SS_scores_Pathogenic.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
envscores <- read.delim("T3SS_scores_Environmental.txt", header=F, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)
rhizscores <- read.delim("T3SS_scores_Rhizosphere.txt", header=F, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

plotbox <- function(gene) {
	set1 <- pathscores2[,gene]
	set2 <- envscores2[,gene]
	set3 <- rhizscores2[,gene]
	pooled <- c(set1,set2,set3)
	png(paste("boxplots/T3SS/", gene, ".png", sep=""))
	plot(x=seq(0.9,1.1,0.2/(length(set1)-1)), y=set1, ylim=c(min(pooled, na.rm=T)-25, max(pooled, na.rm=T)+25), xlim=c(0.5,3.5), main=paste("Score distributions for ", gene, sep=""), xaxt="n", xlab="", ylab="Bitscore", pch=16, cex=1.3, cex.lab=1.3, cex.main=1.3, cex.axis=1.3, frame=F)
	points(x=seq(1.9,2.1,0.2/(length(set2)-1)), y=set2, pch=16)
	points(x=seq(2.9,3.1,0.2/(length(set3)-1)), y=set3, pch=16)
	axis(side=1, at=c(1,2,3), labels=c("Pathogenic", "Environmental", "Rhizosphere"), cex.axis=1.3)
	dev.off()
}

genes <- row.names(pathscores)
genes <- genes[order(genes)]
#for (i in 1:length(genes)) {
#	gene <- genes[i]
#	plotbox(gene)
#}

png("T3SSgenes.avr.png", width=1000, height=500)
par(mar=c(7,5,3,5), cex=1.5)
plot(1,0,xlim=c(0,9), ylim=c(0,1500), col="white", xaxt="n", ylab="Bitscore", xlab="", cex.axis=1.4, cex.lab=1.4)
for (i in 1:10) {
	gene <- genes[i]
	points(x=rep(i-1.2, length(pathscores2[!is.na(pathscores2[,gene]),gene])), y=pathscores2[!is.na(pathscores2[,gene]),gene], pch=16, col="black")
	points(x=rep(i-1, length(rhizscores2[!is.na(rhizscores2[,gene]),gene])), y=rhizscores2[!is.na(rhizscores2[,gene]),gene], pch=16, col="navy")
	points(x=rep(i-0.8, length(envscores2[!is.na(envscores2[,gene]),gene])), y=envscores2[!is.na(envscores2[,gene]),gene], pch=16, col="coral4")
}
axis(side=1, at=(0:9), labels=genes[0:10], las=2, cex.axis=1.4)
dev.off()

png("T3SSgenes.hop1.png", width=1000, height=500)
par(mar=c(7,5,3,5), cex=1.5)
plot(1,0,xlim=c(0,11), ylim=c(0,1500), col="white", xaxt="n", ylab="Bitscore", xlab="", cex.axis=1.4, cex.lab=1.4)
for (i in 11:22) {
	gene <- genes[i]
	points(x=rep(i-11.2, length(pathscores2[!is.na(pathscores2[,gene]),gene])), y=pathscores2[!is.na(pathscores2[,gene]),gene], pch=16, col="black")
	points(x=rep(i-11, length(rhizscores2[!is.na(rhizscores2[,gene]),gene])), y=rhizscores2[!is.na(rhizscores2[,gene]),gene], pch=16, col="navy")
	points(x=rep(i-10.8, length(envscores2[!is.na(envscores2[,gene]),gene])), y=envscores2[!is.na(envscores2[,gene]),gene], pch=16, col="coral4")
}
axis(side=1, at=(0:11), labels=genes[11:22], las=2, cex.axis=1.4)
dev.off()

png("T3SSgenes.hop2.png", width=1000, height=500)
par(mar=c(7,5,3,5), cex=1.5)
plot(1,0,xlim=c(0,13), ylim=c(0,1500), col="white", xaxt="n", ylab="Bitscore", xlab="", cex.axis=1.4, cex.lab=1.4)
for (i in 23:length(genes)) {
	gene <- genes[i]
	points(x=rep(i-23.2, length(pathscores2[!is.na(pathscores2[,gene]),gene])), y=pathscores2[!is.na(pathscores2[,gene]),gene], pch=16, col="black")
	points(x=rep(i-23, length(rhizscores2[!is.na(rhizscores2[,gene]),gene])), y=rhizscores2[!is.na(rhizscores2[,gene]),gene], pch=16, col="navy")
	points(x=rep(i-22.8, length(envscores2[!is.na(envscores2[,gene]),gene])), y=envscores2[!is.na(envscores2[,gene]),gene], pch=16, col="coral4")
}
axis(side=1, at=(0:13), labels=genes[23:length(genes)], las=2, cex.axis=1.4)
dev.off()