invbaselogit <-
function(baselogit){
	    k <- length(baselogit)+1
	    kminus1 <- k-1
	    probs <- rep(0,k)
 	    probs[k] <- 1/(1+sum(exp(baselogit)))
	    for (i in 1:kminus1){
	    	probs[i] <- probs[k]*exp(baselogit[i])
 	    }
	    list(probs=probs)
}

