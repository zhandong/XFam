\name{WPGMSim}
\alias{WPGMSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Winsorized PGM Gibbs Simulator}
\description{Winsorized PGM Gibbs Sampler (both positive and negatvie relationships)}
\usage{
WPGMSim(n, p, R, alpha, Theta, maxit = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size}
  \item{p}{variable size}
  \item{R}{threshold value for truncating}
  \item{alpha}{a px1 vector}
  \item{Theta}{a pxp symmetric matrix (only off diags matter).}
  \item{maxit}{iterations for Gibbs sampler, default to 10000}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{X}{a nxp data matrix}
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
wpgm.sim <- WPGMSim(10, 3, 2, rep(0.5, 3), matrix(-1, 3,3))
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, p, R, alpha, Theta, maxit = 10000) 
{
    X = matrix(rpois(n * p, 1), n, p)
    iter = 1
    while (iter < maxit) {
        for (j in 1:p) {
            num = exp(matrix(1, n, 1) \%*\% t(alpha[j] * c(0:R) - 
                log(factorial(c(0:R)))) + matrix(c(0:R) \%x\% X[, 
                -j] \%*\% Theta[-j, j], n, R + 1))
            Pmat = num/matrix(apply(num, 1, sum), n, R + 1)
            X[, j] = apply(apply(Pmat, 1, mymult) == 1, 2, which) - 
                1
        }
        iter = iter + 1
    }
    return(X)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
