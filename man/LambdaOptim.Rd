\name{LambdaOptim}
\Rdversion{1.1}
\alias{LambdaOptim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Obtains values for the prior specification for lambda
%%  ~~function to do ... ~~
}
\description{In the single \eqn{f} model we may parameterize in terms of the allele frequencies and \eqn{\lambda=\log((f-f_{\min})/(1-f))}  where \eqn{f_{\min}=-p_{\min}/(1-p_{\min})} and \eqn{p_{\min}} is the minimum allele frequency. The prior for \eqn{\lambda} is assumed normal and this function finds the mean and standard deviation of this normal, given two values for \eqn{f}, with associated probabilities.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LambdaOptim(nsim, bvec, f1, f2, p1, p2, init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{the optimization is carried out by simulating from the joint prior on allele frequencies and \eqn{\lambda}, and this argument gives the number of simulations to take from the prior
%%     ~~Describe \code{nsim} here~~
}
  \item{bvec}{vector of length \eqn{k} of prior specification for the HWE Dirichlet prior, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{bvec} here~~
}
  \item{f1}{first quantile for inbreeding coefficient \eqn{f}
%%     ~~Describe \code{f1} here~~
}
  \item{f2}{second quantile for inbreeding coefficient \eqn{f}
%%     ~~Describe \code{f2} here~~
}
  \item{p1}{probability associated with \code{f1}
%%     ~~Describe \code{p1} here~~
}
  \item{p2}{probability associated with \code{f2}
%%     ~~Describe \code{p2} here~~
}
  \item{init}{initial values for \code{lambdamu} and \code{lambdasd}
%%     ~~Describe \code{init} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{lambdamu}{prior mean for \eqn{\lambda}}
\item{lambdasd}{prior standard deviation for \eqn{\lambda}}
%% ...
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with 
\section{Warning }{This function can be unstable and good starting values may be needed. It is also recommended to check the output by simulating from the given prior to see if the empirical quantiles match with those desired; the function \code{SinglefPrior} may be used for this} 

\seealso{\code{HWEImportSamp}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
bvec <- c(1,1,1,1)
init <- c(-3,log(1.1))
lampr <- LambdaOptim(nsim=10000,bvec=bvec,f1=0,f2=0.26,p1=0.5,p2=0.95,init)
}
\keyword{ htest }
