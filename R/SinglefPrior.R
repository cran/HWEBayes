SinglefPrior <-
function(nsim,alpha,lambdamu,lambdasd){
	    p <- rdirichlet(n=nsim,alpha=alpha)
	    lambda <- rnorm(nsim,mean=lambdamu,sd=lambdasd)
	    lgts <- baselogit(p)$baselogit
	    pmin <- apply(p,1,min)
	    fmin <- -pmin/(1-pmin)
	    f <- (exp(lambda)+fmin)/(exp(lambda)+1)
	    list(p=p,f=f,lgts=lgts,lambda=lambda)
}

