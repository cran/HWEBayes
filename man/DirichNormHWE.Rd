\name{DirichNormHWE}
\Rdversion{1.1}
\alias{DirichNormHWE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates the normalizing constant under the HWE model, with a conjugate prior
}
\description{Function to evaluate the normalizing constant given a conjugate Dirichlet prior and the HWE model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DirichNormHWE(nvec, bvec0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{12},..., n_{1k},n_{22} ..., n_{2k},..., n_{kk}}.
%%     ~~Describe \code{nvec} here~~
}
  \item{bvec0}{vector of length \eqn{k} Dirichlet prior parameters, where \eqn{k} is the number of alleles.
}
}
\value{The normalizing constant.
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. 
}
\author{Jon Wakefield (jonno@u.washington).
}

\seealso{\code{DirichNormSat}
}
\examples{
data(DiabRecess)
DirichNormHWE(nvec=DiabRecess,bvec0=rep(1,4))
}
\keyword{ htest }
