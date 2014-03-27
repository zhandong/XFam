\name{WPGM.network}
\alias{WPGM.network}
\title{Poisson network}
\description{Function to compute the poisson network over X}
\usage{
WPGM.network(X, R, nlams, lmin = 0.001, lambda = NULL, parallel = TRUE, ncores = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix (of Poisson)}
  \item{R}{threshold value for truncating}
  \item{nlams}{number of lambdas for regularization path}
  \item{lmin}{minimum lambda value, default to 0.001}
  \item{lambda}{a vector of nlams lambda for whole regularization path, default to NULL}
  \item{parallel}{logical value to indicate if the network build should be run parallelly in multiple threads, default to TRUE}
  \item{ncores}{number of cores to use for parallel execution, default to 4}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of length of the regularization path, each element of the list represent the networks estimated over the regularization path. Each network is encoded in pxp matrix of coefficients. 
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
	\code{\link{WPGM.path.neighborhood}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, nlams, lmin = 0.001, lambda = NULL, parallel = T, 
    ncores = 4) 
{
    if (is.null(lambda)) {
        lmax = lambdaMax(t(X))
        lambda = exp(seq(log(lmax), log(lmin), l = nlams))
    }
    if (nlams != length(lambda)) {
        print("nlams is not equal to lams")
    }
    ghat = c()
    if (nlams > 0) {
        ghat = array(0, dim = c(nrow(X), nrow(X), length(lambda)))
    }
    wrapper <- function(i) {
        fit = WPGM.path.neighborhood(t(X[-i, ]), X[i, ], R, nlams, 
            lambda = lambda, 0)
        fit$beta = as.matrix(fit$Bmat)
        if (i == 1) {
            ghat[i, 2:nrow(X), ] = fit$beta
        }
        else if (i == nrow(X)) {
            ghat[i, 1:(nrow(X) - 1), ] = fit$beta
        }
        else {
            ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), ]
            ghat[i, (i + 1):nrow(X), ] = fit$beta[i:nrow(fit$beta), 
                ]
        }
        return(ghat[i, , ])
    }
    ghat2 = c()
    if (parallel) {
        library(multicore)
        ghat2 = mclapply(1:nrow(X), wrapper, mc.cores = ncores)
    }
    else {
        ghat2 = lapply(1:nrow(X), wrapper)
    }
    for (i in 1:nrow(X)) {
        ghat[i, , ] = ghat2[[i]]
    }
    ghat = lapply(1:nlams, function(r) {
        return(ghat[, , r])
    })
    return(ghat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
