\name{LambdaPriorChoice}
\Rdversion{1.1}
\alias{LambdaPriorChoice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Called by LambaOptim
%%  ~~function to do ... ~~
}
\description{Internal function, should not be needed. It is the function that is minmized when the mean and standard deviation for \eqn{\lambda} are being found.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LambdaPriorChoice(x, nsim, bvec, f1, f2, p1, p2, init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Proposal for \code{lambdamu} and exp(\code{lambdasd}).
%%     ~~Describe \code{x} here~~
}
  \item{nsim}{number of points to simulate from the prior
%%     ~~Describe \code{nsim} here~~
}
  \item{bvec}{\eqn{k} vector of Dirichlet prior parameters.
%%     ~~Describe \code{bvec} here~~
}
  \item{f1}{first quantile for inbreeding coefficient \eqn{f}
%%     ~~Describe \code{f1} here~~
}
  \item{f2}{second quantile for inbreeding coefficient \eqn{f}
%%     ~~Describe \code{f2} here~~
}
  \item{p1}{probability associated with \code{f1}
%%     ~~Describe \code{p1} here~~
}
  \item{p2}{probability associated with \code{f2}
%%     ~~Describe \code{p2} here~~
}
  \item{init}{initial values for \code{lambdamu} and \code{lambdasd}
%%     ~~Describe \code{init} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns the sum of squares that is being minimized.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
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

\seealso{\code{LambdaOptim}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
