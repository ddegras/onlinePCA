snlpca <- function (Q, x, gamma, center, type = c("exact", "nn")) 
{
	if (!missing(center)) 
		x <- x - center
	type <- match.arg(type)	
	if (type == "exact") {
	    Q <- Q + gamma * x %*% crossprod(x, Q)
	    eig <- eigen(crossprod(Q), TRUE)
	    nonzero <- which(eig$values > sqrt(.Machine$double.eps))
	    iS <- eig$vectors[, nonzero] %*% 
	    	(t(eig$vectors[, nonzero])/sqrt(eig$values[nonzero]))
	    Q <- Q %*% iS
	if (length(nonzero) < ncol(Q)) 
        warning(paste("Matrix 'Q' is not full rank. Returning", 
            length(nonzero), "PC."))
	} else if (type == "nn") {
		y <- crossprod(Q, x)	
		Q <- Q + gamma * (tcrossprod(x,y) - Q %*% tcrossprod(y))
	}
	return(Q)
}
