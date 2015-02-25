comm="enteritidis"
patho="gallinarum"

#fitdistribution <- function(comm, patho) {
	
dbs <- read.table(paste(comm, "-", patho, "/results.dbs",sep=""), header=F)
n <- nrow(dbs[dbs$V10!=0,])

dist <- rnorm(n=n, mean=0, sd=sd(dbs[dbs$V12>0.05,10]))
distdens <- density(dist, bw=0.3)

hist(dbs[dbs$V10!=0,10], breaks=500, xlim=c(-7,7), freq=FALSE, ylim=c(0,0.5))
lines(distdens$x, distdens$y)

#}