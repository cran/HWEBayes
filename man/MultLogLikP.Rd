\name{MultLogLikP}
\Rdversion{1.1}
\alias{MultLogLikP}
\title{Evaluates the Multinomial likelihood under the single f model}
\description{Evaluates the Multinomial likelihood under the single \eqn{f} 
model. The normalizing constant is not included. This function is called by 
a number of other functions, and should not be needed.}
\usage{MultLogLikP(p, f, nvec)}
\arguments{
\item{p}{A set of probabilities, \eqn{p_1,...,p_k}, where \eqn{k} is the number of alleles}
\item{f}{The \eqn{f} parameter}
\item{nvec}{vector of genotype frequencies in the order \eqn{n_{11},
    n_{21}, ..., n_{k1}, n_{22}, ..., n_{k2}, ..., n_{kk}}.}
}
\value{
The unnormalized multinomial log-likelihood.
}
\note{
  \code{MultLogLik} also calculates the multinomial likelihood using an
  alternate parameterization.
}
\references{Wakefield, J. (2010). Bayesian methods for examining Hardy-Weinberg
equilibrium. Biometrics; Vol 66:257-65
  
Weir, B.S. (1996). Genetic Data Analysis II. Sunderland MA: Sinauer.
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{SinglefReject},\code{MultLogLik}}
\keyword{ distribution }
