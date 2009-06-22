\name{HWETriangBF2}
\Rdversion{1.1}
\alias{HWETriangBF2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates the Bayes factor in the k=2 allele case with a "triangular" prior under the null
%%  ~~function to do ... ~~
}
\description{Function to evaluate the Bayes factor \eqn{\Pr(n|\mbox{ HWE })/\Pr(n|\mbox{ saturated model })} in the \eqn{k=2} allele case and with a conjugate (Dirichlet) priors under the saturated model and a "triangular" distribution under the null. the latter is the marginal prior distribution under the (1,1,1) Dirichlet prior under the saturated model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HWETriangBF2(nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{21}, n_{22}}.
%%     ~~Describe \code{nvec} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Bayes factor is returned.
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

\seealso{\code{TriangNormHWE}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
HWETriangBF2(nvec=c(88,10,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
