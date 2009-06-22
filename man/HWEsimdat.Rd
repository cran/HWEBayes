\name{HWEsimdat}
\Rdversion{1.1}
\alias{HWEsimdat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate data under the single f model with k alleles.

%%  ~~function to do ... ~~
}
\description{Simulate data under the single \eqn{f} model with \eqn{k} alleles (so \eqn{f=0} gives data under HWE).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HWEsimdat(npop, q, f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{npop}{population size.
%%     ~~Describe \code{npop} here~~
}
  \item{q}{vector of \eqn{k} allele frequencies.
%%     ~~Describe \code{q} here~~
}
  \item{f}{value of inbreeding coefficient
%%     ~~Describe \code{f} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{nvec}{vector of genotype counts, in the order \eqn{n_{11}, n_{21}, n_{22}, ..., n_{k1}, n_{k2}, ..., n_{kk}}.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Wakefield, J. (2009). Bayesian methods for examining Hardy-Weinberg equilibrium. Biometrics.

Weir, B.S. (1996). Genetic Data Analysis II. Sunderland MA: Sinauer.
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
counts <- HWEsimdat(100,q=c(0.1,0.8,.1),f=0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
