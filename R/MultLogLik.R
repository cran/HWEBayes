MultLogLik <-
function(x,nvec,paramch=1){
	   k <- .5*(-1+sqrt(1+8*length(nvec)))
	   pmarg <- rep(0,k)
	   pmarg <- invbaselogit(x[1:k-1])$probs
	   minp <- min(pmarg)
	   minf <- -minp/(1-minp)
	   if (paramch==1) f <- (exp(x[k])-1)/(exp(x[k])+1)
	   if (paramch==2) {f <- (exp(x[k])+minf)/(exp(x[k])+1)}
#	   cat(x,pmarg,f,minp,"\n")
	   const <- 10^80
	   qvec <- rep(0,k+1)
	   if (is.na(minp)) MultLogLik <- -const
	   if (is.na(f)) MultLogLik <- -const
           if ( is.na(minp) == FALSE & is.na(f) == FALSE & f < -minp/(1-minp) ) {
#	      cat("TROUBLE")
	      MultLogLik <- -const
           }
	   if (is.na(minp) == FALSE & is.na(f) == FALSE & f >= -minp/(1-minp) ){
	      MultLogLik <- 0
	      count <- 1
	      for (i in 1:k){
	       for (j in i:k){
	       	   if (j==i){
		      qvec[count] <- pmarg[i]^2 + f*pmarg[i]*(1-pmarg[i])  
#cat(i,j,count,qvec[count],"\n")
		      count <- count+1
		   }
		   if (j>i){
		      qvec[count] <- 2*pmarg[i]*pmarg[j]*(1-f)	     
#cat(i,j,count,qvec[count],"\n")
			 count <- count+1
		   }
	       }
              }
	      for( i in 1:length(nvec)){
	      	   if (qvec[i]<=0) { #cat("qvec PROBLEMS1",qvec,"\n"); 
		      MultLogLik <- -const}
if (qvec[i]>1)  {cat("PROBLEMS2"); MultLogLik <- -const}
	      	   if (qvec[i]>0&&qvec[i]<1) {
		      #cat(i,qvec[i],"\n")
		      MultLogLik <- MultLogLik + nvec[i]*log(qvec[i])}
              }
	   }
	   list(MultLogLik=MultLogLik)
}

