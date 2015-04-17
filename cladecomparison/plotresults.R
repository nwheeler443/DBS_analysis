path <- read.delim("T3SS_scores_Pathogenic.txt", header=F, stringsAsFactors=F, row.names=1)
path1 <- t(path)
env <- read.delim("T3SS_scores_Environmental.txt", header=F, stringsAsFactors=F, row.names=1)
env1 <- t(env)
rhiz <- read.delim("T3SS_scores_Rhizosphere.txt", header=F, stringsAsFactors=F, row.names=1)
rhiz1 <- t(rhiz)

plot(0,0, xlim=c(0,ncol(path1)+1), ylim=c(0,2000), col="white", xaxt="n")
for (effector in 1:length(colnames(path1))) {
	points(x=rep(effector, length(path1[!is.na(path1[,effector]),effector])), y=path1[!is.na(path1[,effector]),effector])
	points(x=rep(effector, length(env1[!is.na(env1[,effector]),effector])), y=env1[!is.na(env1[,effector]),effector], col="red")
	points(x=rep(effector, length(rhiz1[!is.na(rhiz1[,effector]),effector])), y=rhiz1[!is.na(rhiz1[,effector]),effector], col="blue")
	print(colnames(path1[effector]))
	print(path[effector,])
	print(env[effector,])
	print(rhiz[effector,])
}
axis(side=1, at=(1:length(colnames(path1))), labels=colnames(path1), las=2)