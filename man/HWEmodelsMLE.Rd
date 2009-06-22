\name{HWEmodelsMLE}
\Rdversion{1.1}
\alias{HWEmodelsMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates the maximum likelihood estimates of the parameters of various models in the k allele case
%%  ~~function to do ... ~~
}
\description{Function to obtain the MLEs of parameters under the HWE, single \eqn{f} and saturated models. For the single \eqn{f} model numerical maximization is required if \eqn{k>2}, where \eqn{k} is the number of alleles.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HWEmodelsMLE(nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{vector of genotype frequencies in the order \eqn{n_{11}, n_{21}, n_{22}, ..., n_{k1}, n_{k2}, ..., n_{kk}}.
%%     ~~Describe \code{nvec} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{phat}{matrix of \eqn{k \times k} MLEs of genotype frequencies}
\item{qhat}{MLEs of \eqn{k} allele frequencies under the HWE model}
\item{fqhat}{MLEs of \eqn{k} allele frequencies under the single \eqn{f} model}
\item{fsingle}{MLE of single \eqn{f}}
\item{fmaxloglik}{maximized log-likelihood (without the normalizing constant), under the single \eqn{f} model}
\item{fmin}{estimated lower bound of \eqn{f_{\min}} in the single \eqn{f} model. Under the single \eqn{f} model \eqn{f_{\min} < f <1} where \eqn{f_{\min}=-p_{\min}/(1-p_{\min})} and \eqn{p_{\min}} is the minimum of the allele frequencies.}
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
data(DiabRecess)
HWEmodelsMLE(nvec=DiabRecess)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
