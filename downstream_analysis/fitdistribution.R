comm="enteritidis"
patho="gallinarum"

dbs <- read.table(paste(comm, "-", patho, "/results.dbs",sep=""), header=F)
trim <- dbs[dbs$V10!=0,10]
data <- c(trim, -trim)

library(fitdistrplus)

fit <- fitdistr(data, "t")

dist <- rt(length(trim), df=1.28)
