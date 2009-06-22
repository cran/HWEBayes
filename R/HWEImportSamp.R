HWEImportSamp <-
function(nsim,nvec,ischoice,lambdamu,lambdasd,alpha,
 gmu=rep(0,length(alpha)),gsigma=diag(0,nrow=length(alpha),ncol=length(alpha))){
        PrnH1 <- varterm1 <- priorint <- varprior <- 0
	k <- .5*(-1+sqrt(1+8*length(nvec)))
 	kminus1 <- k-1
	if (ischoice==1) {
	         if(gsigma[1,1]==0) stop("HWImportSamp: You need to supply gmu and gsigma")
# we simulate for the baseline logits and lambda
	         phisamp <- mvrnorm(n=nsim,mu=gmu,Sigma=gsigma)	
	      	 for (i in 1:nsim){
      	      	     pval <- invbaselogit(phisamp[i,1:kminus1])$probs 
		     pmin <- min(pval)
		     fmin <- -pmin/(1-pmin)  
      	 	     test <- MultLogLik(x=phisamp[i,],nvec,paramch=2)$MultLogLik
      	 	     likterm <- test+lfactorial(sum(nvec))-sum(lfactorial(nvec))
#
# Log of the determinant of the (k-1)x(k-1) Jacobean, derivs are:
# partial p_1/partial phi1,...partial p_1/partial phi_{k-1}
# ........
# partial p_{k-1}/partial phi1,...partial p_{k-1}/partial phi_{k-1}
#
		     jac <- matrix(0,nrow=kminus1,ncol=kminus1)
	   	     denom <- (1+sum(exp(phisamp[i,1:kminus1])))^2
	   	     for (j1 in 1:kminus1){
	       		    for (j2 in 1:kminus1){
	       	    	    	if(j1==j2) jac[j1,j1] <- exp(phisamp[i,j1])*
                                         (1+sum(exp(phisamp[i,1:kminus1]))-exp(phisamp[i,j1]))/denom
	   	   		if(j1>j2) jac[j1,j2] <- -exp(phisamp[i,j1]+phisamp[i,j2])/denom
		   		if(j1<j2) jac[j1,j2] <- -exp(phisamp[i,j1]+phisamp[i,j2])/denom
	       	            }
	             }
	             ljack <- log(det(jac))
		     logjac1 <- ljack + phisamp[i,k] + log(1-fmin) - log((1+exp(phisamp[i,k]))^2)
      	 	     prterm1 <- log(ddirichlet(pval, alpha=alpha)) + logjac1
		     f <- (exp(phisamp[i,k])+fmin)/(exp(phisamp[i,k])+1)
	      	     lambda <- log((f-fmin)/(1-f))
	      	     lambdav <- lambdasd^2
	      	     logdens <- -.5*log(2*pi*lambdasd*lambdasd) - .5*(lambda-lambdamu)^2/(lambdasd*lambdasd)
	      	     logjac <- log(1-fmin) - log(f-fmin) - log(1-f)
	      	     prterm2 <- logdens + logjac
      	 	     gterm <- dmvnorm(phisamp[i,1:k],mean=gmu,sigma=gsigma,log=TRUE) 
      	 	     expterm <- exp(likterm+prterm1+prterm2-gterm)
#		     expprior <- exp(prterm1+prterm2-gterm)
      	 	     PrnH1 <- PrnH1 + expterm
#		     priorint <- priorint + expprior
      	 	     varterm1 <- varterm1 + expterm^2
#		     varprior <- varprior + expprior^2
		     if (i/1000 - round(i/1000) == 0) cat("Samples = ",i,"\n")
     	          }
#		  priorint <- priorint/nsim
#     		  varprior <- (varprior/nsim - priorint^2)/nsim
#     		  cat("nsim prior constant (se) 95% interval = ",nsim,priorint,"(",sqrt(varprior),")",priorint-1.96*sqrt(varprior),priorint+1.96*sqrt(varprior),"\n")
    }
     if (ischoice==2){
     	for (i in 1:nsim){
	      pval <- rdirichlet(1,alpha=alpha)
	      minp <- min(pval)
	      fmin <- -minp/(1-minp)
	      lambdaval <- rnorm(1,mean=lambdamu,sd=lambdasd)
	      test <- MultLogLik(x=c(baselogit(pval)$baselogit,lambdaval),nvec,paramch=2)$MultLogLik
      	      likterm <- test+lfactorial(sum(nvec))-sum(lfactorial(nvec))
      	      expterm <- exp(likterm)
      	      PrnH1 <- PrnH1 + expterm
      	      varterm1 <- varterm1 + expterm^2
	 }
     }
     PrnH1 <- PrnH1/nsim
     varest <- (varterm1/nsim - PrnH1^2)/nsim
     cat("nsim norm constant (se) 95% interval = ",nsim,PrnH1,"(",sqrt(varest),")",PrnH1-1.96*sqrt(varest),PrnH1+1.96*sqrt(varest),"\n")
     list(PrnH1=PrnH1,varest=varest)
}

