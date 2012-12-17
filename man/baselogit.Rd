\name{baselogit}
\Rdversion{1.1}
\alias{baselogit}
\title{Calculates a set of baseline logits from a set of probabilities
}
\description{Calculates a set of \eqn{k-1} baseline logits \eqn{\log(p_1/p_k), ..., \log(p_{k-1}/p_k)}, from a set of probabilities \eqn{p_1,...,p_k}.
}
\usage{
baselogit(probs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probs}{A set of probabilities, \eqn{p_1,...,p_k}, where \eqn{k} is the number of alleles.
}
}
\details{This function is used by a number of other functions in the package, for example, to provide a parameterization for maximization and for importance sampling in the single \eqn{f} model.
}
\value{
\item{baselogit}{Returns the set of \eqn{k-1} baseline logits, where \eqn{k} is the number of alleles.}
}
\references{Wakefield, J. (2010). Bayesian methods for examining Hardy-Weinberg
equilibrium. Biometrics; Vol 66:257-65}
\author{Jon Wakefield (jonno@u.washington).
}
\seealso{\code{invbaselogit}
}
\examples{
baselogit(probs=c(0.5,0.4,0.1))
}
\keyword{ misc }
