SinglefReject <-
function(nsim,bvec,lambdamu,lambdasd,nvec){
    k <- .5*(-1+sqrt(1+8*length(nvec)))
    MLEres <- HWEmodelsMLE(nvec)
    maxLL <- MLEres$fmaxloglik
    psamp <- matrix(0,nrow=nsim,ncol=k)
    fsamp <- NULL
    naccept <- count <- PrnH1 <- varterm1 <- 0
    while (naccept < nsim){
      count <- count+1
      samples <- SinglefPrior(nsim=1,alpha=bvec,lambdamu=lambdamu,
                              lambdasd=lambdasd)
      LL <- MultLogLik(c(samples$lgt,samples$lambda),nvec,paramch=2)$MultLogLik
      likterm <- LL+lfactorial(sum(nvec))-sum(lfactorial(nvec))
      expterm <- exp(likterm)
      PrnH1 <- PrnH1 + expterm
      varterm1 <- varterm1 + expterm^2
      if (LL > maxLL) cat("Maximization is messed up\n")
      u <- runif(1)
      if (log(u) < LL-maxLL){
         naccept <- naccept+1
         psamp[naccept,1:k] <- samples$p
         fsamp[naccept] <- samples$f
         if (floor(naccept/100)==ceiling(naccept/100))cat("Number of accepted points = ",naccept,"\n")
      }
   }
   PrnH1 <- PrnH1/count
   varest <- (varterm1/count - PrnH1^2)/count
   cat("nsim norm constant (se) 95% interval = ",nsim,PrnH1,"(",sqrt(varest),")",PrnH1-1.96*sqrt(varest),PrnH1+1.96*sqrt(varest),"\n")	
   accrate <- nsim/count
   list(psamp=psamp,fsamp=fsamp,accrate=accrate,PrnH1=PrnH1,varest=varest)
}

