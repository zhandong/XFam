\name{WPGM.neighborhood}
\alias{WPGM.neighborhood}
\title{WPGM neighborhood}
\description{WPGM neighborhood selection problem (on one lambda)
}
\usage{
WPGM.neighborhood(X, Y, R, lam, startb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a nxp data matrix}
  \item{Y}{nx1 vector of responses (Poisson?)}
  \item{R}{threshold value for truncating}
  \item{lam}{numeric lambda value (regularization parameter)}
  \item{startb}{default to 0, otherwise a starting vector for beta}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of:
  \item{alpha}{intercept}
  \item{beta}{vector of p coefficients}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y, R, lam, startb = 0) 
{
    n = nrow(X)
    p = ncol(X)
    thr = 1e-08
    maxit = 1e+06
    Xt = cbind(t(t(rep(1, n))), X)
    if (sum(startb) == 0) {
        bhat = matrix(rnorm(p + 1) * 0.01, p + 1, 1)
    }
    else {
        bhat = startb
    }
    step = 0.1
    ind = 1
    iter = 1
    while (thr < ind & iter < maxit) {
        oldb = bhat
        t = 1
        grad = wpgmGrad(Xt, Y, R, oldb)
        oldobj = wpgmObj(Xt, Y, R, oldb)
        tmp = oldb - t * grad
        bhat[1] = tmp[1]
        bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - lam * 
            t, max, 0)
        while (wpgmObj(Xt, Y, R, bhat) > oldobj - t(grad) \%*\% 
            (oldb - bhat) + sum((oldb - bhat)^2)/(2 * t)) {
            t = t * step
            tmp = oldb - t * grad
            bhat[1] = tmp[1]
            bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                lam * t, max, 0)
        }
        iter = iter + 1
        ind = sum((oldb - bhat)^2)/sum(oldb^2)
    }
    return(list(alpha = bhat[1], beta = bhat[-1]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
