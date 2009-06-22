\name{DirichSampSat}
\Rdversion{1.1}
\alias{DirichSampSat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate samples from a Dirichlet prior or posterior under the saturated model
%%  ~~function to do ... ~~
}
\description{Function to simulate samples from the satuated Dirichlet model. Can be used for samples from the prior or the (conjugate) Dirichlet posterior, both in the \eqn{k} allele case. Samples are generated for the genotype frequencies in the order \eqn{p_{11}, p_{12},..., p_{1k},p_{22} ..., p_{2k},..., p_{kk}}, the allele frequencies, and the fixation indices.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{DirichSampSat(nvec, bvec, nsim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{12},..., n_{1k},n_{22} ..., n_{2k},..., n_{kk}}.
%%     ~~Describe \code{nvec} here~~
}
  \item{bvec}{vector of length \eqn{k(k+1)/2} Dirichlet prior parameters, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{bvec} here~~
}
  \item{nsim}{number of samples to simulate from the prior/posterior.
%%     ~~Describe \code{nsim} here~~
}
}
\details{Uses the \code{rdirichlet} function from the \code{MCMCpack} library.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{pvec}{matrix of size \code{nsim} \eqn{\times k(k+1)/2} containing samples for the genotype frequencies, in the order \eqn{p_{11}, p_{12},..., p_{1k},p_{22} ..., p_{2k},..., p_{kk}}.}
\item{pmat}{matrix of size \code{nsim} \eqn{\times k(k+1)/2 \times k(k+1)/2} containing samples for the genotype probabilities.}
\item{pmarg}{matrix of size \code{nsim} \eqn{\times k} containing samples for the allele frequencies, in the order \eqn{p_{1},..., p_k}.}
\item{fixind}{matrix of size \code{nsim} \eqn{\times k(k+1)/2 \times k(k+1)/2} 
containing samples for the fixation indices, contained in the lower diagonal, 
i.e. \code{fixind[,i,j]} for \code{[i>j]}.}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{DirichSampHWE}, \code{DirichNormSat}, \code{DirichNormHWE}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# First sample from the prior
PriorSampSat <- DirichSampSat(nvec=rep(0,10),bvec=rep(1,10),nsim=1000)
par(mfrow=c(1,2))
hist(PriorSampSat$pvec[,1],xlab="p1",main="")
hist(PriorSampSat$fixind[,2,1],xlab="f21",main="")
# Now sample from the posterior
data(DiabRecess)
PostSampSat <- DirichSampSat(nvec=DiabRecess,bvec=rep(1,10),nsim=1000)
par(mfrow=c(1,2))
hist(PostSampSat$pvec[,1],xlab="p1",main="")
hist(PostSampSat$fixind[,2,1],xlab="f21",main="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution}

