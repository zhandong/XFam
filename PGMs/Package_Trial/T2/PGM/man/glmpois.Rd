\name{glmpois}
\alias{glmpois}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Poisson based neighborhood selection}
\description{Poisson based neighborhood selection with X (pxn) on a fixed regularization path. }
\usage{
glmpois(X, lambda, parallel = F, nCpus = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{lambda}{regularization parameter, could be a single numerical value or a vector of numeric values (for the whole regularization path)}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to FALSE}
  \item{nCpus}{number of (maximum) cores to be used for parallel execution, default to 4}
}
\details{
This function will depends on the \code{glmnet} function.  The neighborhood selection method is based on Meinshausen and Buhlmann neighborhood selection methods proposed for Gaussian graphical models.
}
\value{
  \item{ghat}{If a specific lambda is given, ghat is a pxp matrix of coefficients. 
			If the lambda for the whole regularization path is input, a 3D (pxpx length of regularization path) matrix is returned, where \code{ghat[,,i]} is the coefficient matrix of the p variables for the \code{i}-lambda
		}
}
\references{
N. Meinshausen and P. Buhlmann, 2006, High-dimensional graphs and variable selection with the lasso, \emph{The Annals of Statistics}, \bold{vol. 34}, no. 3, pp. 1436--1462.

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link[glmnet]{glmnet}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, lambda, parallel = F, nCpus = 4) 
{
    if (length(lambda) > 1) {
        ghat = array(0, dim = c(nrow(X), nrow(X), length(lambda)))
        if (parallel) {
            wrapper <- function(i) {
                fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                  lambda = lambda)
                fit$beta = as.matrix(fit$beta)
                if (ncol(fit$beta) < length(lambda)) {
                  tmp = matrix(0, nrow = nrow(fit$beta), ncol = length(lambda))
                  tmp[, 1:ncol(fit$beta)] = fit$beta
                  tmp[, ncol(fit$beta):length(lambda)] = fit$beta[, 
                    ncol(fit$beta)]
                  fit$beta = tmp
                }
                if (i == 1) {
                  ghat[i, 2:nrow(X), ] = fit$beta
                }
                else if (i == nrow(X)) {
                  ghat[i, 1:(nrow(X) - 1), ] = fit$beta
                }
                else {
                  ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), 
                    ]
                  ghat[i, (i + 1):nrow(X), ] = fit$beta[i:nrow(fit$beta), 
                    ]
                }
                return(ghat[i, , ])
            }
            library(multicore)
            ghat2 = mclapply(1:nrow(X), wrapper)
            for (i in 1:nrow(X)) {
                ghat[i, , ] = ghat2[[i]]
            }
            return(ghat)
        }
        if (parallel == F) {
            wrapper <- function(i) {
                fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                  lambda = lambda)
                fit$beta = as.matrix(fit$beta)
                if (ncol(fit$beta) < length(lambda)) {
                  tmp = matrix(0, nrow = nrow(fit$beta), ncol = length(lambda))
                  tmp[, 1:ncol(fit$beta)] = fit$beta
                  tmp[, ncol(fit$beta):length(lambda)] = fit$beta[, 
                    ncol(fit$beta)]
                  fit$beta = tmp
                }
                if (i == 1) {
                  ghat[i, 2:nrow(X), ] = fit$beta
                }
                else if (i == nrow(X)) {
                  ghat[i, 1:(nrow(X) - 1), ] = fit$bet
                }
                else {
                  ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), 
                    ]
                  ghat[i, (i + 1):nrow(X), ] = fit$beta[i:nrow(fit$beta), 
                    ]
                }
                return(ghat[i, , ])
            }
            ghat2 = lapply(1:nrow(X), wrapper)
            for (i in 1:nrow(X)) {
                ghat[i, , ] = ghat2[[i]]
            }
            return(ghat)
        }
    }
    if (length(lambda) == 1) {
        ghat = matrix(0, nrow = nrow(X), ncol = nrow(X))
        if (parallel) {
            library(snowfall)
            sfInit(cpus = nCpus)
            sfExport("X", local = T)
            sfExport("ghat", local = T)
            sfLibrary(glmnet)
            wrapper <- function(i) {
                fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                  lambda = lambda)
                fit$beta = as.numeric(fit$beta)
                if (i == 1) {
                  ghat[i, 2:nrow(X)] = fit$beta
                }
                else if (i == nrow(X)) {
                  ghat[i, 1:(nrow(X) - 1)] = fit$beta
                }
                else {
                  ghat[i, 1:(i - 1)] = fit$beta[1:(i - 1)]
                  ghat[i, (i + 1):nrow(X)] = c(fit$beta[i:length(fit$beta)])
                }
                return(ghat[i, ])
            }
            sfExport("wrapper")
            ghat = sfSapply(1:nrow(X), wrapper)
            sfStop()
            return(ghat)
        }
        for (i in 1:nrow(X)) {
            fit = glmnet(t(X[-i, ]), X[i, ], family = "poisson", 
                lambda = lambda)
            fit$beta = as.numeric(fit$beta)
            if (i == 1) {
                ghat[i, 2:nrow(X)] = fit$beta
            }
            else if (i == nrow(X)) {
                ghat[i, 1:(nrow(X) - 1)] = fit$beta
            }
            else {
                ghat[i, 1:(i - 1)] = fit$beta[1:(i - 1)]
                ghat[i, (i + 1):nrow(X)] = c(fit$beta[i:length(fit$beta)])
            }
        }
        return(ghat)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
