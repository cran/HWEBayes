\name{HWEDirichBF2}
\Rdversion{1.1}
\alias{HWEDirichBF2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates the Bayes factor in the k=2 allele case under conjugate priors
}
\description{Function to evaluate the Bayes factor \eqn{\Pr(n|\mbox{ HWE })/\Pr(n|\mbox{ saturated model })} in the \eqn{k=2} allele case and with conjugate (Dirichlet) priors under HWE and saturated models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HWEDirichBF2(nvec, bvec0, bvec1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{21}, n_{22}}.
%%     ~~Describe \code{nvec} here~~
}
  \item{bvec0}{vector of length \eqn{k=2} Dirichlet prior parameters for the prior under the null, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{bvec0} here~~
}
  \item{bvec1}{vector of length \eqn{k(k+1)/2=3} Dirichlet prior parameters for the prior under the saturated model, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{bvec1} here~~
}
}
\value{Bayes factor is returned.
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{DirichNormHWE}, \code{DirichNormSat}, \code{DirichSampHWE}, \code{DirichSampSat}, \code{HWETriangBF2}, \code{TriangNormHWE}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
HWEDirichBF2(nvec=c(88,10,2),bvec0=c(1,1),bvec1=c(1,1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
