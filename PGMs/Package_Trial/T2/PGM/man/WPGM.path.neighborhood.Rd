\name{WPGM.path.neighborhood}
\alias{WPGM.path.neighborhood}
\title{WPGM neighborhood over a regularization path}
\description{WPGM neighborhood selection problem over a grid of lambdas
}
\usage{
WPGM.path.neighborhood(X, Y, R, nlams, lmin = 0.01, lambda = NULL, startb = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a nxp data matrix}
  \item{Y}{nx1 vector of responses (Poisson?)}
  \item{R}{threshold value for truncating}
  \item{nlams}{number of lambdas for regularization path (set nlams=1 to return form one value)}
  \item{lmin}{minimum lambda value, default to 0.01}
  \item{lambda}{a vector of nlams lambda, default to NULL}
  \item{startb}{default to 0, otherwise a starting vector for beta}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of:
  \item{alphas}{1 x nlams vector of intercepts}
  \item{Bmat}{p x nlams sparse matrix of coefficients}
  \item{lambda}{the lambda values for regularization path}
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
function (X, Y, R, nlams, lmin = 0.01, lambda = NULL, startb = 0) 
{
    n = nrow(X)
    p = ncol(X)
    if (is.null(lambda)) {
        lmax = lambdaMax(t(X))
        lambda = exp(seq(log(lmax), log(lmin), l = nlams))
    }
    if (nlams == 1 & is.null(lambda)) {
        lambda = lmax
    }
    thr = 1e-08
    maxit = 1e+06
    Xt = cbind(t(t(rep(1, n))), X)
    if (sum(startb) == 0) {
        bhat = matrix(rnorm(p + 1)/p, p + 1, 1)
    }
    else {
        bhat = startb
    }
    alphas = 0
    Bmat = matrix(0, p, nlams)
    step = 0.1
    for (i in 1:nlams) {
        ind = 1
        iter = 1
        while (thr < ind & iter < maxit) {
            oldb = bhat
            t = 1
            grad = wpgmGrad(Xt, Y, R, oldb)
            oldobj = wpgmObj(Xt, Y, R, oldb)
            tmp = oldb - t * grad
            bhat[1] = tmp[1]
            bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                lambda[i] * t, max, 0)
            newobj = wpgmObj(Xt, Y, R, bhat)
            while (newobj > 9999999 | is.na(newobj) | is.na(newobj)) {
                t = t/p
                tmp = oldb - t * grad
                bhat[1] = tmp[1]
                bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                  lambda[i] * t, max, 0)
                newobj = wpgmObj(Xt, Y, R, bhat)
            }
            while (newobj > oldobj - t(grad) \%*\% (oldb - bhat) + 
                sum((oldb - bhat)^2)/(2 * t)) {
                t = t * step
                tmp = oldb - t * grad
                bhat[1] = tmp[1]
                bhat[-1] = sign(tmp[-1]) * sapply(abs(tmp[-1]) - 
                  lambda[i] * t, max, 0)
                newobj = wpgmObj(Xt, Y, R, bhat)
            }
            iter = iter + 1
            ind = sum((oldb - bhat)^2)
        }
        alphas[i] = bhat[1]
        Bmat[, i] = bhat[-1]
    }
    return(list(alpha = alphas, Bmat = Bmat, lambda = lambda))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
