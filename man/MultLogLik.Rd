\name{MultLogLik}
\Rdversion{1.1}
\alias{MultLogLik}
\title{Evaluates the Multinomial likelihood under the single f model}
\description{Evaluates the Multinomial likelihood under the single \eqn{f} 
model. The normalizing constant is not included. This function is called by 
a number of other functions, and should not be needed.}
\usage{MultLogLik(x, nvec, paramch = 1)}
\arguments{
\item{x}{a set of \eqn{k-1} baseline logits, where \eqn{k} is the number of 
alleles), and a transformed version of \eqn{f}. Hence a vector of length 
\eqn{k}. The transformation adopted depends on the value of \code{paramch}.}
\item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{21}, n_{22}, ..., n_{k1}, n_{k2}, ..., n_{kk}}.}
\item{paramch}{a variable that if =1 assumes \eqn{f} is on the range (-1,+1) before transformation, and if =2 assumes on the range \eqn{(f_{\min},+1}).}
}
\details{}
\value{
\item{MultLoglik}{The value of the  (unnormalized) multinomial log-likelihood.}
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.

Weir, B.S. (1996). Genetic Data Analysis II. Sunderland MA: Sinauer.
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{SinglefReject}}
\examples{}
\keyword{ distribution }
