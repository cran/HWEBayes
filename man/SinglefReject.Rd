\name{SinglefReject}
\Rdversion{1.1}
\alias{SinglefReject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Samples from the posterior for the single f model
}
\description{Function to generate samples from the posterior for allele frequencies and \eqn{f}, under the single \eqn{f} model.
Samples are generated using a rejection algorithm that simulates from the prior.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SinglefReject(nsim, bvec, lambdamu, lambdasd, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{number of samples to generate from the prior.
%%     ~~Describe \code{nsim} here~~
}
  \item{bvec}{vector of size \eqn{k} that is the specification for the Dirichlet prior on the allele frequencies.
%%     ~~Describe \code{bvec} here~~
}
  \item{lambdamu}{prior mean for \eqn{\lambda}.
%%     ~~Describe \code{lambdamu} here~~
}
  \item{lambdasd}{prior standard deviation for \eqn{\lambda}.
%%     ~~Describe \code{lambdasd} here~~
}
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{21}, n_{22}, ..., n_{k1}, n_{k2}, ..., n_{kk}}.
%%     ~~Describe \code{nvec} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{psamp}{samples for \eqn{k} allele frequencies.}
\item{fsamp}{samples for inbreeding coefficient \eqn{f}.}
\item{accrate}{acceptance rate of the rejection algorithm.}
\item{PrnH1}{estimate of normalizing constant (which may be used in Bayes factor calculations). Calculated by averaging the likelihood over the sampled points.}
\item{varest}{estimated variance of the estimate of the normalizing constant.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Wakefield, J. (2010). Bayesian methods for examining Hardy-Weinberg
equilibrium. Biometrics; Vol 66:257-65
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(DiabRecess)
postsampf1 <- SinglefReject(nsim=100,bvec=rep(1,4),lambdamu=-2.95,
   lambdasd=1.07,nvec=DiabRecess)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
