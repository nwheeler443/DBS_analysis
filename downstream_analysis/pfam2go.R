pfam2go <- function(comm, patho){

#comm <- "typhimurium"
#patho <- "typhi"

	pfams <- read.delim("~/Documents/Thesis/Other functional mapping/pfam2go.extended.fix", header=F, stringsAsFactors=F)
	
	dbs <- read.table(paste(comm, "-", patho, "/results.dbs", sep=""))

	dbs[,3] <- as.factor(sub("\\.\\d+", "", dbs[,3]))
	length <- nrow(dbs)

	hash <- function( keys ) {
		result <- new.env( hash = TRUE, parent = emptyenv(), size = length( keys ) )
		for( key in keys ) {
			result[[ key ]] <- NA
		}
		return( result )
	}

	lines <- c()

	pval_hash <- hash(unique(pfams$V3))

	for(i in unique(pfams$V3)) {
		this_term <- pfams[grepl(i, pfams$V3),]
		labelled_list <- cbind(dbs, rep(0, nrow(dbs)))
		for(j in this_term[,1]) {
			labelled_list[grep(j, labelled_list[,3]),13]=1
	#print(j)
		}
		labelled_list <- labelled_list[order(labelled_list[,10], decreasing=TRUE),]
		hits <- sum(labelled_list[, 13])
		n_lab <- sum(labelled_list[,13])
		n_unlab <- length - n_lab
		
		max_p <- 0
		max_m <- 0
		vals <- c()
		if (hits > 0) {
			for (j in 1:length){
				hit <- sum(labelled_list[1:j, 13])
		#calculate hypergeometric p-value for observing this many or more at this position
				p_val_up <- phyper(hit - 1, n_lab, n_unlab, j, lower.tail=FALSE)
				p_val_down <- phyper(n_lab - hit - 1, n_lab, n_unlab, length - j, lower.tail=FALSE)
				log_p_up <- -log(p_val_up, base=10)
				log_p_down <- -log(p_val_down, base=10)
				if(labelled_list[j,10] >= 0){
	#print("success")
					if(log_p_up > max_p){
						max_p <- log_p_up
						max_p_stat <- labelled_list[j,10]
						max_p_g <- hit
						max_p_g_up <- j
					}
				}
				else if(labelled_list[j,10] < 0){
					if(log_p_down > max_m){
						max_m <- log_p_down
						max_m_stat <- labelled_list[j,10]
						max_m_g <- n_lab - hit
						max_m_g_down <- length - j
					}
				}
				vals <- rbind(vals, c(j, log_p_up, log_p_down))
		#}
			}
			
		#calculate correlation with phenotype
			up_cor <- cor(labelled_list[,10],vals[,2], method="spearman")
			down_cor <- cor(-labelled_list[,10], vals[,3], method="spearman")
			
			pval_hash[[i]] <- vals
			
			lines <- rbind(lines, c(as.character(this_term[1,2]),as.character(this_term[1,3]), n_lab ,max_p,  up_cor, max_p_stat, max_p_g,max_p_g_up, max_m, down_cor, max_m_stat,max_m_g,max_m_g_down))
		}
	}

  lines[,4] <- -log(p.adjust(10^(-as.numeric(lines[,4])), method="BH"), base=10)
  lines[,9] <- -log(p.adjust(10^(-as.numeric(lines[,9])), method="BH"), base=10)
  
  lines <- lines[order(apply(lines,1, function(row) max(as.numeric(row[4]), as.numeric(row[9]))), decreasing=TRUE), ]
  colnames(lines) <- c("pathway","pathname","path_genes", "up_p","up_cor", "up_max_stat","path_genes_up", "genes_up", "down_p","down_cor","down_max_stat","path_genes_down", "genes_down")
  write.table(lines, file=paste(comm, "-", patho, "/goanalysis.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)

}