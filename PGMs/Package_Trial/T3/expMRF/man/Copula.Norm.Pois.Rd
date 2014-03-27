\name{Copula.Norm.Pois}
\alias{Copula.Norm.Pois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Copula transform a matrix from normal to Poisson}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Copula.Norm.Pois(X, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a nxp data matrix of Gaussians}
  \item{lambda}{the Poisson mean for the transformation}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Y}{a nxp Copula transformed data matrix}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
		X <- matrix(rnorm(20), nrow=5, ncol=4)
		transX <- Copula.Norm.Pois(X, lambda=1)

## The function is currently defined as
function (X, lambda) 
{
    n = nrow(X)
    p = ncol(X)
    val = 0
    dcuts = NULL
    cnt = 0
    while (val < max(0.9999, 1 - 2/(n * p))) {
        val = ppois(cnt, lambda)
        cnt = cnt + 1
        dcuts = c(dcuts, val)
    }
    Y = matrix(0, n, p)
    oldval = min(X)
    for (i in 1:length(dcuts)) {
        val = quantile(X, dcuts[i])
        Y[which(X < val & X >= oldval)] = i - 1
        oldval = val
    }
    Y[X == max(X)] = max(Y) + 1
    return(Y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
