baselogit <-
function(probs){
	  k <- length(probs)
	  kminus1 <- k-1; baselogit <- rep(0,kminus1)
	  for (i in 1:kminus1){
	      baselogit[i] <- log(probs[i]/probs[k])}
	  list(baselogit=baselogit)
}

