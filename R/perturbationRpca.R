perturbationRpca <- function (d, Q, x, n, ff, center) 
{
    if (missing(ff)) 
        ff <- (1 / n) else if (ff <= 0 || ff >= 1) 
			stop("Argument 'ff' must be in (0,1)")
	p <- length(x) 
    if (length(d) != p)
    	stop("'d' and 'x' of incompatible dimensions")
    if (ncol(Q) != p)
    	stop("'d' and 'Q' of incompatible dimensions")
    if (!missing(center)) 
    	x <- x - center
    newd <- (1 - ff) * d
	newd2 <- newd^2
    newz <- sqrt((1 + ff) * ff) * crossprod(Q, x)
    newz2 <- newz^2
	num <- tcrossprod(newz)  
	den <- matrix(newd + newz2, p, p, byrow = TRUE) - 
		matrix(newz2 + newd2, p, p)
	V <- num / den
	diag(V) <- 1
	Q <- Q %*% V
    sigma2 <- .colSums(Q^2, p, p)
    newQ <- Q / rep(sqrt(sigma2), each = p)
	newd <- (newd + newz2) * sigma2
    list(values = newd, vectors = newQ)
}
