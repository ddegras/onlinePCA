sgapca <- function(Q, x, gamma, center, type = c("exact", "nn"))  
{
	if (!missing(center)) 
		x <- x - center
	type <- match.arg(type)
	if (type == "exact") {
	    y <- as.vector(crossprod(Q, x))
	    Q <- Q + tcrossprod(x, gamma * y)
	    Q <- qr.Q(qr(Q))
	} else if (type == "nn") {
	    y <- rep(crossprod(Q, x), each = nrow(Q))
		M <- diag(ncol(Q))
		M[upper.tri(M)] <- 2
		M[lower.tri(M)] <- 0
		Q <- Q + gamma * y * (x - {y * Q} %*% M)
	}
	return(Q)
}
