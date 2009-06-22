\name{HWEImportSamp}
\Rdversion{1.1}
\alias{HWEImportSamp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Importance sampling to calculate the normalizing constant under the single f model
}
\description{Importance sampling to calculate the normalizing constant under the single \eqn{f} model. Two proposals are available, either sampling from the prior or sampling from a normal distribution whose mean vector and variance-covariance matrix must be specified. The latter may be taken from an MCMC analysis using, for example, WinBUGS. In all cases the likelihood is multinomial and the prior is Dirichlet on the allele frequencies, and normal on \eqn{\lambda} where \eqn{\lambda=\log(f-f_{\min})/(1-f))}. See Weir (1996) for a description of HWE and different models/parameterizations.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HWEImportSamp(nsim, nvec, ischoice, lambdamu, lambdasd, alpha,
gmu = rep(0, length(alpha)), gsigma = diag(0, nrow = length(alpha), 
ncol = length(alpha)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{the number of points to sample to calculate the estimate.
%%     ~~Describe \code{nsim} here~~
}
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{12},..., n_{1k},n_{22} ..., n_{2k},..., n_{kk}}.
%%     ~~Describe \code{nvec} here~~
}
  \item{ischoice}{choice of importance sampling proposal, =1 gives a normal distribution with mean and variance that must be specified (as \code{gmu} and \code{gsigma}) and =2 is from the prior.
%%     ~~Describe \code{ischoice} here~~
}
  \item{lambdamu}{the mean of the prior for \eqn{\lambda}.
%%     ~~Describe \code{lambdamu} here~~
}
  \item{lambdasd}{the variance of the prior for \eqn{\lambda}.
%%     ~~Describe \code{lambdasd} here~~
}
  \item{alpha}{the vector of \eqn{k} parameters for the Dirichlet prior on the allele frequencies.
%%     ~~Describe \code{alpha} here~~
}
  \item{gmu}{the mean of the importance sampling proposal, of length \eqn{k}, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{gmu} here~~
}
  \item{gsigma}{the variance of the importance sampling proposal, a matrix of dimension \eqn{k \times k}, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{gsigma} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{PrnH1}{the estimate of the normalizing constant}
\item{varest}{the variance of the estimate of the normalizing constant}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.

Weir, B.S. (1996). Genetic Data Analysis II. Sunderland MA: Sinauer.
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{As always with importance sampling the procedure can be very unstable, particularly for large \eqn{k}. Hence rerunning the function with different simulation sample sizes, and different \code{gmu} and \code{gsigma} is recommended}
\seealso{\code{LambdaOptim}, \code{DirichNormSat}, \code{DirichNormHWE}, \code{TriangNormHWE}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
alpha <- c(1,1,1,1) # prior on allele frequencies
# gmu and gsigma were obtained from a WinBUGS run
gmu <- c(-0.4633092,0.3391625,0.3397936,-3.5438008)
gsigma <- matrix(c(
0.07937341,0.02819656,0.02766583,0.04607996,
0.02819656,0.07091320,0.04023827,0.01657028,
0.02766583,0.04023827,0.07042278,0.01752266,
0.04607996,0.01657028,0.01752266,0.57273683),nrow=4,ncol=4)
data(DiabRecess)
HWEImportSamp(nsim=5000,nvec=DiabRecess,ischoice=1,lambdamu=-2.95,
   lambdasd=1.07,alpha=alpha,gmu,gsigma)
HWEImportSamp(nsim=5000,nvec=DiabRecess,ischoice=2,lambdamu=-2.95,
   lambdasd=1.07,alpha=alpha)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
