LambdaPriorChoice <-
function(x,nsim,bvec,f1,f2,p1,p2,init){
	   p <- rdirichlet(nsim,alpha=bvec)	
	   pmin <- apply(p,1,min)
	   fmin <- -pmin/(1-pmin)	
	   lambdamu <- x[1]
	   lambdasd <- exp(x[2])
	   lambda <- rnorm(nsim,lambdamu,lambdasd)
           f <- (exp(lambda)+fmin)/(exp(lambda)+1)
	   (quantile(f,probs=p1)-f1)^2+(quantile(f,probs=p2)-f2)^2
}

