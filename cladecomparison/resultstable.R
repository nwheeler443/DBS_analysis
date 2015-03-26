pathenv <- read.delim("path-env.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("pathogenscores.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
envscores <- read.delim("environmentalscores.txt", header=F, row.names=1, stringsAsFactors=F)
envscores2 <- t(envscores)

kstable <- data.frame()
for (orth in 1:nrow(pathenv)) {
	if (!is.na(match(pathenv[orth,1], colnames(pathscores2)))) {
		set1 <- pathscores2[,pathenv[orth,1]]
		set2 <- envscores2[,pathenv[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		if(ksval$p.value<0.005) {
			png(paste("boxplots/", pathenv[orth,1], ".png", sep=""))
			boxplot(set1,set2)
			dev.off()
		}
		ksrow <- data.frame(gene.name=pathenv[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, p.adj=padj)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(pathscores))

write.table(kstable, file="path-env.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05,1], file="path-env.siggenes", quote=F, col.names=F, row.names=F, sep="\t")

pathrhiz <- read.delim("path-rhiz.dbs/orthlist.dbs", header=F, stringsAsFactors=F)

pathscores <- read.delim("pathogenscores.txt", header=F, row.names=1, stringsAsFactors=F)
pathscores2 <- t(pathscores)
rhizscores <- read.delim("rhizospherescores.txt", header=F, row.names=1, stringsAsFactors=F)
rhizscores2 <- t(rhizscores)

kstable <- data.frame()
for (orth in 1:nrow(pathrhiz)) {
	if (!is.na(match(pathrhiz[orth,1], colnames(pathscores2))) & !is.na(match(pathrhiz[orth,2], colnames(rhizscores2)))) {
		set1 <- pathscores2[,pathrhiz[orth,1]]
		set2 <- rhizscores2[,pathrhiz[orth,2]]
		ksval = ks.test(set1, set2, exact=FALSE)
		if(ksval$p.value<0.005) {
			png(paste("boxplots/rhiz.", pathrhiz[orth,1], ".png", sep=""))
			boxplot(set1,set2)
			dev.off()
		}
		ksrow <- data.frame(gene.name=pathrhiz[orth,1], statistic=as.numeric(ksval$statistic), p.val=ksval$p.value, p.adj=padj)
		kstable <- rbind(kstable, ksrow)
	}
}

kstable$p.adj = p.adjust(kstable$p.val, method="BH", n=nrow(pathscores))

write.table(kstable, file="path-rhiz.ksvalues", quote=F, row.names=F, col.names=T, sep="\t")
write.table(kstable[kstable$p.adj<0.05,1], file="path-rhiz.siggenes", quote=F, col.names=F, row.names=F, sep="\t")
