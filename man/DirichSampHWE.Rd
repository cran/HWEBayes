\name{DirichSampHWE}
\Rdversion{1.1}
\alias{DirichSampHWE}
\title{Simulate samples from a Dirichlet prior or posterior under HWE 
}
\description{Function to simulate samples from the HWE Dirichlet model. Can be used for samples from the prior or the (conjugate) Dirichlet posterior, both in the \eqn{k} allele case. Samples are generated for the allele frequencies in the order \eqn{p_{1},p_{2},...,p_{k}}.
}
\usage{DirichSampHWE(nvec, bvec0, nsim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{12},..., n_{1k},n_{22} ..., n_{2k},..., n_{kk}}.
%%     ~~Describe \code{nvec} here~~
}
  \item{bvec0}{vector of length \eqn{k} Dirichlet prior parameters, where \eqn{k} is the number of alleles.
}
  \item{nsim}{number of samples to simulate from the prior/posterior.
}
}
\details{Uses the \code{rdirichlet} function from the \code{MCMCpack} library.
}
\value{
\item{pvec}{matrix of size \code{nsim} \eqn{\times k} containing samples for the genotype frequencies, in the order \eqn{p_{1}, p_{12},..., p_{k}}.}
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.
}
\author{Jon Wakefield (jonno@u.washington).
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{DirichSampSat}, \code{DirichNormSat}, \code{DirichNormHWE}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# First sample from the prior
PriorSampHWE <- DirichSampHWE(nvec=rep(0,10),bvec0=rep(1,4),nsim=1000)
par(mfrow=c(1,1))
hist(PriorSampHWE$pvec[,1],xlab="p1",main="")
# Now sample from the posterior
data(DiabRecess)
PostSampHWE <- DirichSampHWE(nvec=DiabRecess,bvec0=rep(1,4),nsim=1000)
par(mfrow=c(1,1))
hist(PostSampHWE$pvec[,1],xlab="p1",main="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
