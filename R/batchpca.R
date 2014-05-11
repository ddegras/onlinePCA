batchpca <- function(C, k)
{
	if (missing(k))
		k <- ncol(C)
	if (k <= ncol(C) / 10) {
		res <- eigs_sym(C, k) 
	} else { res <- eigen(C, TRUE)
		res$values <- res$values[seq_len(k)]
		res$vectors <- res$vectors[,seq_len(k)] 
	}	
	list(values = res$values, vectors = res$vectors)
}
