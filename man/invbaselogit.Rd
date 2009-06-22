\name{invbaselogit}
\Rdversion{1.1}
\alias{invbaselogit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converts a set of k-1 baseline logits into a set of probabilities
%%  ~~function to do ... ~~
}
\description{Converts a set of \eqn{k-1} baseline logits \eqn{\log (p_1/p_k),...,\log (p_{k-1}/p_k)} into a set of probabilities \eqn{p_1,...,p_k}, where \eqn{k} is the number of alleles. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
invbaselogit(baselogit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baselogit}{A set of \eqn{k-1} baseline logits, where \eqn{k} is the number of alleles.
%%     ~~Describe \code{baselogit} here~~
}
}
\details{This is used by a number of other functions in the package, for example, to provide a parameterization for maximization and for importance sampling in the single \eqn{f} model.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{invbaselogit }{the probability vector corresponding to the baseline logit}
%%  \item{comp2 }{Description of 'comp2'}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{baselogit}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
invbaselogit(baselogit=c(0,0))
}
\keyword{ misc }
