updateCovariance <- function (C, x, n, xbar, ff, byrow = TRUE) 
{
	if (missing(n) && missing(ff)) 
		stop("At least one of the arguments 'n' and 'ff' must be specified")
	if (!is.matrix(x)) {
        x <- as.matrix(x)
        byrow <- FALSE
    }
    dimc <- dim(C)
	if (dimc[1] != dimc[2]) 
		stop("'C' must be a square matrix")
	p <- dimc[1]
    dimx <- dim(x)
    k <- ifelse(byrow, dimx[1], dimx[2]) 
    pp <- ifelse(byrow, dimx[2], dimx[1]) 
    if (p != pp) 
        stop(paste0("'C' and 'x' of incompatible dimensions.\n",
        	"Check these arguments and 'byrow'"))
	if (!missing(xbar) && length(xbar) != p) 
       stop(paste0("'x' and 'xbar' of incompatible dimensions.\n",
        	"Check these arguments and 'byrow'"))

    meanx <- if (byrow) {
    	.colMeans(x, k, p) 
    	} else .rowMeans(x, p, k)
	if (missing(ff)) 
		ff <- 1 / (n + k - 1)	
    ffm <- 1 / (1 + 1 / ff)
   
    if (missing(xbar)) {
    	newxbar = Dxbar = (k * ffm) * meanx 
    	} else {
    	newxbar <- (1 - k * ffm) * xbar + (k * ffm) * meanx 
    	Dxbar <- newxbar - xbar
    }
	a1 <- 1 - k * ff
	a2 <- 1 - (k - 1) * ff
	x <- x - matrix(newxbar, dimx[1], dimx[2], byrow)

    if (k == 1L) 
        return(a1 * C + a2 * tcrossprod(Dxbar))
	if (byrow)
		return(a1 * C + a2 * tcrossprod(Dxbar) + 
			ff * crossprod(x))
	return(a1 * C + a2 * tcrossprod(Dxbar) + 
		ff * tcrossprod(x))
}