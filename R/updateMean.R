updateMean <- function(xbar, x, n, ff, byrow = TRUE)
{
	if (missing(n) && missing(ff)) 
		stop("At least one of the arguments 'n' and 'ff' must be specified")
	if (!is.matrix(x)) {
		x <- as.matrix(x)
		byrow <- FALSE}
	dimx <- dim(x)
    k <- ifelse(byrow, dimx[1], dimx[2]) 
    p <- ifelse(byrow, dimx[2], dimx[1]) 
	if (length(xbar) != p) 
       stop(paste0("'x' and 'xbar' of incompatible dimensions.\n",
        	"Check these arguments and 'byrow'"))
	if (missing(ff)) 
		ff <- k / (n + k)
	mean_x <- if (byrow) {
		.colMeans(x, k, p) } else .rowMeans(x, p, k)		
	return((1 - ff) * xbar + ff * mean_x)
}