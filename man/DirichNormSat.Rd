\name{DirichNormSat}
\Rdversion{1.1}
\alias{DirichNormSat}
\title{Evaluates the normalizing constant (as used in the denominator of a Bayes factor) for a conjugate prior
}
\description{Function to evaluate the normalizing constant given a conjugate Dirichlet prior and a saturated model.
}
\usage{
DirichNormSat(nvec, bvec)
}
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{12},..., n_{1k},n_{22} ..., n_{2k},..., n_{kk}}.
}
  \item{bvec}{vector of length \eqn{k(k+1)/2} Dirichlet prior parameters, where \eqn{k} is the number of alleles.
}
}
\value{The normalizing constant.
}
\references{Wakefield, J. (2010). Bayesian methods for examining Hardy-Weinberg
equilibrium. Biometrics; Vol 66:257-65
}
\author{Jon Wakefield (jonno@u.washington.edu)}

\seealso{\code{DirichSampSat}, \code{DirichSampHWE}, \code{DirichNormSat}, \code{DirichNormHWE}, \code{HWEDirichBF2}
}
\examples{
data(DiabRecess)
DirichNormSat(nvec=DiabRecess,bvec=rep(1,10))
}
\keyword{ htest }
