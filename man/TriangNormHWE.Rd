\name{TriangNormHWE}
\Rdversion{1.1}
\alias{TriangNormHWE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates the normalizing constant under the HWE model, for the "triangular" prior distribution
}
\description{Function to evaluate the normalizing constant given a "triangular" prior and the HWE model, in the \eqn{k=2} allele case. This prior results from marginalizing the conjugate Dirichlet prior with parameters (1,1,1) on the genotype frequencies under the alternative. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TriangNormHWE(nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{21}, n_{22}}.
%%     ~~Describe \code{nvec} here~~
}
}
\value{Normalizing constant is returned.
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{DirichNormHWE}, \code{DirichNormSat}, \code{DirichSampHWE}, \code{DirichSampSat}, \code{HWEDirichBF2}, \code{TriangNormHWE}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nvec <- c(88,10,2)
TriangNormHWE(nvec)
}
\keyword{ htest }
