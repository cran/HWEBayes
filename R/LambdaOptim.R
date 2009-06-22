LambdaOptim <-
function(nsim,bvec,f1,f2,p1,p2,init){
	    opt <- optim(par=init,fn=LambdaPriorChoice,nsim=nsim,bvec=bvec,
	    	f1=f1,f2=f2,p1=p1,p2=p2)
            lambdamu <- opt$par[1]
	    lambdasd = exp(opt$par[2])
	    cat("lambda mu and lambda sd = ",lambdamu,lambdasd,"\n")
	    list(lambdamu=lambdamu,lambdasd=lambdasd)
}

