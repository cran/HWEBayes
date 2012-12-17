\name{SinglefPrior}
\Rdversion{1.1}
\alias{SinglefPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Samples from the single f prior.
}
\description{Function to sample from the single \eqn{f} prior, that is the Dirichlet and normal on \eqn{\lambda}, where \eqn{\lambda=\log((f-f_{\min})/(1-f))}. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SinglefPrior(nsim, alpha, lambdamu, lambdasd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{number of simulations from prior
%%     ~~Describe \code{nsim} here~~
}
  \item{alpha}{vector of \eqn{k} parameters for the Dirichlet prior on the \eqn{k} allele frequencies.
%%     ~~Describe \code{alpha} here~~
}
  \item{lambdamu}{mean of the normal prior on \eqn{\lambda}.
%%     ~~Describe \code{lambdamu} here~~
}
  \item{lambdasd}{standard deviation of the normal prior on \eqn{\lambda}.
%%     ~~Describe \code{lambdasd} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{p}{sample for vector of \eqn{k} allele frequencies}
\item{f}{sample of \eqn{f} parameters}
\item{lgts}{samples for logits of baseline logits}
\item{lambda}{samples for \eqn{\lambda}}
%% ...
}
\references{Wakefield, J. (2010). Bayesian methods for examining Hardy-Weinberg
equilibrium. Biometrics; Vol 66:257-65
%% ~put references to the literature/web site here ~
}
\author{Jon Wakefield (jonno@u.washington.edu)
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{SinglefReject}, \code{HWEsimdat}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{SinglefSamp <- SinglefPrior(nsim=1000,alpha=c(1,1,1,1),
  lambdamu=-2.95,lambdasd=1.07)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
