pathenv <- read.delim("path-env.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("pathogenscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
envscores <- read.delim("environmentalscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)

kstable <- data.frame()
for (orth in 1:nrow(pathenv)) {
	if (!is.na(match(pathenv[orth,1], colnames(pathscores2)))) {
		set1 <- pathscores2[,pathenv[orth,1]]
		set2 <- envscores2[,pathenv[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		DBS = median(set1, na.rm=T)-median(set2, na.rm=T)
		if(ksval$statistic>0.9&abs(DBS)>50) {
			png(paste("boxplots/NAs/", pathenv[orth,1], ".png", sep=""))
			boxplot(set1,set2, names=c("pathogenic", "environmental"))
			dev.off()
		}
		ksrow <- data.frame(gene.name=pathenv[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, DBS=DBS)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(pathscores))
kstable$z.score = (kstable$DBS-0)/sd(kstable$DBS)
kstable$DBS.p.value = 2*pnorm(-abs(kstable$z.score))
kstable$DBS.p.adj = p.adjust(kstable$DBS.p.val, method="BH", n=nrow(pathscores))

write.table(kstable, file="path-env.NAs.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05&kstable$DBS.p.adj<0.05,1], file="path-env.NAs.siggenes", quote=F, col.names=F, row.names=F, sep="\t")

pathrhiz <- read.delim("path-rhiz.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("pathogenscores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
rhizscores <- read.delim("rhizospherescores.NAs.txt", header=F, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

kstable <- data.frame()
for (orth in 1:nrow(pathrhiz)) {
	if (!is.na(match(pathrhiz[orth,1], colnames(pathscores2))) & !is.na(match(pathrhiz[orth,2], colnames(rhizscores2)))) {
		set1 <- pathscores2[,pathrhiz[orth,1]]
		set2 <- rhizscores2[,pathrhiz[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		DBS = median(set1, na.rm=T)-median(set2, na.rm=T)
		if(ksval$statistic>0.9&abs(DBS)>50) {
			png(paste("boxplots/NAs/", pathrhiz[orth,1], ".png", sep=""))
			boxplot(set1,set2, names=c("pathogenic", "rhizosphere"))
			dev.off()
		}
		ksrow <- data.frame(gene.name=pathrhiz[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, DBS=DBS)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(pathscores))
kstable$z.score = (kstable$DBS-0)/sd(kstable$DBS)
kstable$DBS.p.value = 2*pnorm(-abs(kstable$z.score))
kstable$DBS.p.adj = p.adjust(kstable$DBS.p.val, method="BH", n=nrow(pathscores))

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
		if(ksval$statistic>0.9&abs(DBS)>50) {
			png(paste("boxplots/NAs/", envrhiz[orth,1], ".png", sep=""))
			boxplot(set1,set2, names=c("environmental", "rhizosphere"))
			dev.off()
		}
		ksrow <- data.frame(gene.name=envrhiz[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, DBS=DBS)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(envscores))
kstable$z.score = (kstable$DBS-0)/sd(kstable$DBS)
kstable$DBS.p.value = 2*pnorm(-abs(kstable$z.score))
kstable$DBS.p.adj = p.adjust(kstable$DBS.p.val, method="BH", n=nrow(pathscores))

write.table(kstable, file="env-rhiz.NAs.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05&kstable$DBS.p.adj<0.05,1], file="env-rhiz.NAs.siggenes", quote=F, col.names=F, row.names=F, sep="\t")

