pathenv <- read.delim("path-env.dbs/orthlist.dbs", header=F, stringsAsFactors=F)
pathrhiz <- read.delim("path-rhiz.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("pathogenscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
envscores <- read.delim("environmentalscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)
rhizscores <- read.delim("rhizospherescores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

path_env <- read.delim("path-env.NAs.genes", header=F, stringsAsFactors=F)
path_rhiz <- read.delim("path-rhiz.NAs.genes", header=F, stringsAsFactors=F)

kstable <- data.frame()
for (orth in 1:nrow(pathenv)) {
	if (!is.na(match(pathenv[orth,1], colnames(pathscores2)))) {
		set1 <- pathscores2[,pathenv[orth,1]]
		set2 <- envscores2[,pathenv[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		DBS = median(set1, na.rm=T)-median(set2, na.rm=T)
#		if(ksval$statistic>0.9&abs(DBS)>50) {
#			png(paste("boxplots/NAs/", pathenv[orth,1], ".png", sep=""))
#			boxplot(set1,set2, names=c("pathogenic", "environmental"))
#			dev.off()
#		}
		ksrow <- data.frame(gene.name=pathenv[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, DBS=DBS)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(pathscores))
dbstrim <- kstable$DBS[kstable$DBS<(mean(kstable$DBS)+4*sd(kstable$DBS)) & kstable$DBS>(mean(kstable$DBS)-4*sd(kstable$DBS))]
lengths <- length(dbstrim)
dbstrim <- dbstrim[dbstrim<(mean(dbstrim)+4*sd(dbstrim)) & dbstrim>(mean(dbstrim)-4*sd(dbstrim))]
while (min(lengths) != length(dbstrim)) {
	lengths <- append(lengths, length(dbstrim))
	dbstrim <- dbstrim[dbstrim<(mean(dbstrim)+4*sd(dbstrim)) & dbstrim>(mean(dbstrim)-4*sd(dbstrim))]
}
kstable$z.score = (kstable$DBS-0)/sd(dbstrim)
kstable$DBS.p.value = 2*pnorm(-abs(kstable$z.score))
kstable$DBS.p.adj = p.adjust(kstable$DBS.p.val, method="BH", n=nrow(pathscores))
kstable <- kstable[,c(1,2,3,5,4,6,7,8)]

dist <- rnorm(n=500000, mean=0, sd=sd(dbstrim))
distdens <- density(dist, bw=0.3)

png("path-env.png", width=500, height=500)
hist(kstable$DBS, breaks=500, freq=F, xlim=c(-50,50), main="Score distribution for Pathogenic-Environmental comparison", xlab="Delta-bitscore", ylab="Density")
lines(distdens, col="coral4", lwd=2)
dev.off()

write.table(kstable, file="path-env.NAs.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05&kstable$DBS.p.adj<0.05,1], file="path-env.NAs.siggenes", quote=F, col.names=F, row.names=F, sep="\t")

kstable <- data.frame()
for (orth in 1:nrow(pathrhiz)) {
	if (!is.na(match(pathrhiz[orth,1], colnames(pathscores2))) & !is.na(match(pathrhiz[orth,2], colnames(rhizscores2)))) {
		set1 <- pathscores2[,pathrhiz[orth,1]]
		set2 <- rhizscores2[,pathrhiz[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		DBS = median(set1, na.rm=T)-median(set2, na.rm=T)
#		if(ksval$statistic>0.9&abs(DBS)>50) {
#			png(paste("boxplots/NAs/", pathrhiz[orth,1], ".png", sep=""))
#			boxplot(set1,set2, names=c("pathogenic", "rhizosphere"))
#			dev.off()
#		}
		ksrow <- data.frame(gene.name=pathrhiz[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, DBS=DBS)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(pathscores))
dbstrim <- kstable$DBS[kstable$DBS<(mean(kstable$DBS)+4*sd(kstable$DBS)) & kstable$DBS>(mean(kstable$DBS)-4*sd(kstable$DBS))]
lengths <- length(dbstrim)
dbstrim <- dbstrim[dbstrim<(mean(dbstrim)+4*sd(dbstrim)) & dbstrim>(mean(dbstrim)-4*sd(dbstrim))]
while (min(lengths) != length(dbstrim)) {
	lengths <- append(lengths, length(dbstrim))
	dbstrim <- dbstrim[dbstrim<(mean(dbstrim)+4*sd(dbstrim)) & dbstrim>(mean(dbstrim)-4*sd(dbstrim))]
}
kstable$z.score = (kstable$DBS-0)/sd(dbstrim)
kstable$DBS.p.value = 2*pnorm(-abs(kstable$z.score))
kstable$DBS.p.adj = p.adjust(kstable$DBS.p.val, method="BH", n=nrow(pathscores))
kstable <- kstable[,c(1,2,3,5,4,6,7,8)]

dist <- rnorm(n=500000, mean=0, sd=sd(dbstrim))
distdens <- density(dist, bw=0.3)

png("path-rhiz.png", width=500, height=500)
hist(kstable$DBS, breaks=500, freq=F, xlim=c(-50,50), main="Score distribution for Pathogenic-Rhizosphere comparison", xlab="Delta-bitscore", ylab="Density")
lines(distdens, col="coral4", lwd=2)
dev.off()

write.table(kstable, file="path-rhiz.NAs.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05&kstable$DBS.p.adj<0.05,1], file="path-rhiz.NAs.siggenes", quote=F, col.names=F, row.names=F, sep="\t")

envrhiz <- read.table("env-rhiz.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

envscores <- read.delim("environmentalscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)
rhizscores <- read.delim("rhizospherescores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

kstable <- data.frame()
for (orth in 1:nrow(envrhiz)) {
	if (!is.na(match(envrhiz[orth,1], colnames(envscores2))) & !is.na(match(envrhiz[orth,2], colnames(rhizscores2)))) {
		set1 <- envscores2[,envrhiz[orth,1]]
		set2 <- rhizscores2[,envrhiz[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		DBS = median(set1, na.rm=T)-median(set2, na.rm=T)
#		if(ksval$statistic>0.9&abs(DBS)>50) {
#			png(paste("boxplots/NAs/", envrhiz[orth,1], ".png", sep=""))
#			boxplot(set1,set2, names=c("environmental", "rhizosphere"))
#			dev.off()
#		}
		ksrow <- data.frame(gene.name=envrhiz[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, DBS=DBS)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(envscores))
dbstrim <- kstable$DBS[kstable$DBS<(mean(kstable$DBS)+4*sd(kstable$DBS)) & kstable$DBS>(mean(kstable$DBS)-4*sd(kstable$DBS))]
lengths <- length(dbstrim)
dbstrim <- dbstrim[dbstrim<(mean(dbstrim)+4*sd(dbstrim)) & dbstrim>(mean(dbstrim)-4*sd(dbstrim))]
while (min(lengths) != length(dbstrim)) {
	lengths <- append(lengths, length(dbstrim))
	dbstrim <- dbstrim[dbstrim<(mean(dbstrim)+4*sd(dbstrim)) & dbstrim>(mean(dbstrim)-4*sd(dbstrim))]
}
kstable$z.score = (kstable$DBS-0)/sd(dbstrim)
kstable$DBS.p.value = 2*pnorm(-abs(kstable$z.score))
kstable$DBS.p.adj = p.adjust(kstable$DBS.p.val, method="BH", n=nrow(pathscores))
kstable <- kstable[,c(1,2,3,5,4,6,7,8)]

dist <- rnorm(n=500000, mean=0, sd=sd(dbstrim))
distdens <- density(dist, bw=0.3)

png("rhiz-env.png", width=500, height=500)
hist(kstable$DBS, breaks=500, freq=F, xlim=c(-50,50), main="Score distribution for Rhizosphere-Environmental comparison", xlab="Delta-bitscore", ylab="Density")
lines(distdens, col="coral4", lwd=2)
dev.off()

write.table(kstable, file="env-rhiz.NAs.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05&kstable$DBS.p.adj<0.05,1], file="env-rhiz.NAs.siggenes", quote=F, col.names=F, row.names=F, sep="\t")

