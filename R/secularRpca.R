secularRpca <- function (d, Q, x, n, ff, center, tol = 1e-10, reortho = FALSE) 
{
 
	if (missing(ff)) 
		ff <- 1 / n else if (ff <= 0 || ff >= 1) 
			stop("Argument 'ff' must be in (0,1)")
    if (!missing(center)) 
    	x <- x - center
	p <- length(d)
	if (p != ncol(Q)) 
		stop("Arguments 'd' and 'Q' of incompatible dimensions")	
    d <- (1 - ff) * d
    z <- sqrt((1 + ff) * ff) * crossprod(Q, x)
	ix <- order(d)
	Q <- Q[,ix]
	d <- d[ix]
	z <- z[ix]
	eps <- max(tol, .Machine$double.eps)
	active <- seq_len(p)
		
	# 1.1. Initial deflation: null components of the update vector
	zzero <- which(abs(z) < eps)
	if (length(zzero)) active <- setdiff(active, zzero)
	
	# 1.2. Initial deflation: multiple eigenvalues 
	# of the current eigendecomposition
	ind <- which(diff(d[active]) < eps)
	if (length(ind)) {
		ixequal <- active[sort(unique(c(ind,ind+1L)))]
		dequal <- d[ixequal]	# multiple eigenvalues (w/ repetition)
		u <- unique(dequal)		# multiple eigenvalues (w/o repetition)
		nmult <- length(u)		# number of distinct multiple eigenvalues 
		breaks <- c(-Inf, (u[-nmult] + u[-1]) / 2, Inf)
		group <- split(ixequal, findInterval(dequal, breaks))
		# Each component of the list 'group' contains the indexes 
		# of a multiple eigenvalue 
		
		for (i in seq_along(group))
		{
			g <- group[[i]]
			mult <- length(g)	# multiplicity of eigenvalue
			Q1 <- Q[,g]
			z1 <- z[g]
			sigma <- sqrt(sum(z1^2))
			v <- c(sigma + z1[1], z1[-1])
			v <- v / sqrt(sum(v^2))
			H <- diag(length(g)) - 2 * tcrossprod(v)	# elementary reflector
			Q[,g] <- tcrossprod(Q1, H)
			z[g] <- c(-sigma, rep(0,mult-1))
			active <- setdiff(active,g[-1])
			rm(g, Q1, z1, sigma, v, H)
		}	
	}
		# At this point, the active pairs (d,z) are such that 
		# all d's are unique and all z's are nonzero

	# 2. Computation of eigenvalues via secular equation 
	pact <- length(active)
	dact <- d[active]
	z2act <- z[active]^2
	bounds <- c(dact, dact[pact] + sum(z2act))	
	amp <- diff(bounds)	
	f <- function(lambda) sum(z2act / {dact - lambda}) + 1
	
	solver <- function(i)
	{
		delta <- amp[i] / 100
		repeat {
			lb <- bounds[i] + delta
			ub <- bounds[i+1] - delta
			flb <- f(lb)
			fub <- f(ub) 
			test <- flb * fub

			if (test < 0) {
				return(uniroot(f, c(lb, ub), f.lower = flb, f.upper = fub, 
				tol = tol)$root)
			} else if (test == 0) {
				return(ifelse(flb == 0, lb, ub))
			} else if (flb > 0 && lb - bounds[i] <= tol) {
				return((bounds[i] + lb) / 2)
			} else if (fub < 0 && bounds[i+1] - ub <= tol) {
				return((ub + bounds[i+1]) / 2)
			} else delta <- delta / 10			
		}		
	}

	roots <- sapply(seq_len(pact), solver)	
	d[active] <- roots
	 
	# 3. Associated eigenvectors
	if (reortho)
		{ num <- matrix(roots - rep(dact, each = pact), pact, pact)
		den <- matrix(dact - rep(dact, each = pact), pact, pact)	
		den[den==0] <- 1
		ztilde <- sqrt(apply(num/den, 2, prod))
		eigvecs <- matrix(ztilde, pact, pact) / 
			matrix(dact - rep(roots, each = pact), pact, pact)
		norms <- sqrt(.colSums(eigvecs^2, pact, pact) )
		eigvecs <- eigvecs / rep(norms, each = pact)
	} else {
		eigvecs <- matrix(z[active], pact, pact) / 
			matrix(dact - rep(roots, each = pact), pact, pact)  
		norms <- sqrt(.colSums(eigvecs^2, pact, pact) )
		eigvecs <- eigvecs / rep(norms, each = pact)
	}

	Q[,active] <- Q[,active] %*% eigvecs
	return(list(values = d[p:1], vectors = Q[,p:1]))
}
