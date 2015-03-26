stackedbar <- function(comm, path) {

	pathdata<-read.table(paste(comm, "-", path, "/pathloss.genes.unique.txt", sep=""), header = T, sep = "\t")
	pathdata <- pathdata[nrow(pathdata):1,]
	
	pathdata$totalgenes <- pathdata[,1]*2+(pathdata[,4]+pathdata[,5])*2
	pathdata[,7:10] <- pathdata[,2:5]/pathdata$totalgenes						# losses and uniques as fractions of total
pathdata[,11:12] <- (pathdata[,1]-pathdata[,2:3]+pathdata[,4:5])/pathdata[,6]			# functional proportion of genome
	#pathdata <- pathdata[order(pathdata[,13]),]	

	pathdata <- pathdata[pathdata[,12]!="Inf",]
	pathdataT<-matrix(c(pathdata[,10], pathdata[,7], pathdata[,11], pathdata[,12], pathdata[,8], pathdata[,9]),nrow=nrow(pathdata),ncol=6,dimnames=list(c(as.character(rownames(pathdata))),c("unique to pathogen", "loss in commensal", "commensal", "pathogen", "loss in pathogen", "unique to commensal")))
	
	simpleCap <- function(x) {
		s <- strsplit(x, " ")[[1]]
		paste(toupper(substring(s, 1,1)), substring(s, 2),		# this capitalizes all words - maybe delete one of the 1s to fix this?
			  sep="", collapse=" ")
	}
	
	png(file=paste(comm, "-", path, "/pathlosses.png", sep=""), width=1000, height=500)
	op<-par(las=2, xpd=T, mar=c(0,20,5,20))				# mfrow for mulit-panelled plot
	bp <- barplot(t(pathdataT) , col=c(rgb(30, 30, 100, max=255), rgb(70, 95, 190, max=255), rgb(130, 179, 255, max=255), rgb(205, 82, 82, max=255), rgb(145, 42, 42, max=255), rgb(100, 0, 0, max=255)),main=paste("Genes in ", comm, " vs ", path, sep=""),ylab="Loss (%)", horiz=TRUE, cex.lab=0.2, xaxt="n", cex.names=1.3, cex.main=1.5)
	#bp <- barplot(t(pathdataT) , col=c(rgb(30, 30, 100, max=255), rgb(70, 95, 190, max=255), rgb(130, 179, 255, max=255), rgb(205, 82, 82, max=255), rgb(145, 42, 42, max=255), rgb(100, 0, 0, max=255)),main=paste("Genes in Enteritidis vs Gallinarum", sep=""),ylab="Loss (%)", horiz=TRUE, cex.lab=0.2, xaxt="n", cex.names=1.3, cex.main=1.5)
	
	text(x=c(rep(.45, nrow(pathdataT)), rep(.55, nrow(pathdataT))), y=bp, labels=c(pathdata[,6]/2-pathdata[,5]-pathdata[,2],pathdata[,6]/2-pathdata[,4]-pathdata[,3]), cex=1.1)
	
	legend(x=1.025, y=nrow(pathdata)/2+3, cex=1.25, c(paste("loss in ", comm, "/gain in ", path, sep=""), paste("loss of function in ", comm, sep=""), paste("functional in ", comm, sep=""), paste("functional in ", path, sep=""), paste("loss of function in ", path, sep=""), paste("loss in ", path, "/gain in ", comm, sep="")), fill=c(rgb(30, 30, 100, max=255), rgb(70, 95, 190, max=255), rgb(130, 179, 255, max=255), rgb(205, 82, 82, max=255), rgb(145, 42, 42, max=255), rgb(100, 0, 0, max=255)))
#abline(v=0.5)
#text(x=c(.25,.75), y=nrow(pathdataT)+5, labels=c("Control", "Bad bug"), font=2, cex=1.3)
#text(x=c(.25,.75), y=nrow(pathdataT)+5, labels=c("Enteritidis", "Gallinarum"), font=2, cex=1.3)
text(x=c(.25,.75), y=nrow(pathdataT)+5, labels=c(comm, path), font=2, cex=1.3)
	
	dev.off()
	
	pdf(file=paste(comm, "-", path, "/pathlosses.pdf", sep=""), width=12, height=7)
	op<-par(las=2, xpd=T, mar=c(0,20,5,2))
	bp <- barplot(t(pathdataT) , col=c(rgb(30, 30, 100, max=255), rgb(70, 95, 190, max=255), rgb(130, 179, 255, max=255), rgb(205, 82, 82, max=255), rgb(145, 42, 42, max=255), rgb(100, 0, 0, max=255)),main=paste("Functional changes in ", comm, " vs ", path, sep=""),horiz=TRUE, cex.lab=1.2, cex.axis=1.5, cex.names=1.2, cex.main=1.5, xaxt="n")
	#bp <- barplot(t(pathdataT) , col=c("blue4","cornflowerblue", "brown1", "brown4"),main="Functional changes in Enteritidis vs Gallinarum",horiz=TRUE, cex.lab=1.2, cex.axis=1.5, cex.names=1.5, cex.main=2.5, xaxt="n")
	text(x=c(rep(25, nrow(pathdataT)), rep(75, nrow(pathdataT))), y=bp, labels=c(pathdata[,1],pathdata[,2]), cex=1.1,pos=2)
	text(x=c(rep(33, nrow(pathdataT)), rep(83, nrow(pathdataT))), y=bp, labels=c(paste("(", pathdata[,3], ")", sep=""), paste("(", pathdata[,4], ")", sep="")), cex=1.1,pos=2)
	text(x=c(25,75), y=nrow(pathdataT)+5, labels=c(comm, path), font=4, cex=1.2)
	#text(x=c(25,75), y=nrow(pathdataT)+5, labels=c("Enteritidis", "Gallinarum"), font=2, cex=1.2)
	#text(x=50, y=-4, labels="Gene content (%)", cex=1.2)
	par(xpd=F)
#abline(v=50, lty=2)
	dev.off()
	
}