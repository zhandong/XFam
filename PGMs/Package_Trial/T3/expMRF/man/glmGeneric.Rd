\name{glmGeneric}
\alias{glmGeneric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic Function for local log-linear graphical model.
}
\description{
	(Generic) Function to implement the  local Log-Linear Graphical Model based on pair-wise markov network with the efficient and parallel algorithm.
}
\usage{
	glmGeneric(X, Y = NULL, link, lambda, parallel = FALSE, nCpus = 4, standardize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
	a pxn data matrix
}
  \item{Y}{
	a qxn data matrix or NULL, defualt to NULL. If it's a data matrix, the column number (n) should be the same as X. Y should be of different distribution family as in X.
}
  \item{link}{
	link family to specify the distribution family of the response. If Y is NULL, link indicates the family of X. If Y is not NULL, link indicates the family of Y.
}
  \item{lambda}{
	lambda vector used for regularization
}
  \item{parallel}{
	logical value to indicate if the process should be run in parallel on multiple threads, default to FALSE.
}
  \item{nCpus}{
	number of (maximum) cores to use for parallel execution, default to 4.
}
  \item{standardize}{
	Logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is
standardize=TRUE. This parameter passed to glmnet for the parameter of the same name.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Y = NULL, link, lambda, parallel = F, nCpus = 4, 
    standardize = TRUE) 
{
    if (is.null(Y)) {
        Z <- X
        p <- nrow(Z)
        q <- 0
    }
    if (!is.null(Y)) {
        if (ncol(X) == ncol(Y)) {
            Z <- rbind(X, Y)
            p = nrow(X)
            q = nrow(Y)
        }
    }
    if (length(lambda) > 1) {
        ghat = array(0, dim = c(nrow(Z), nrow(Z), length(lambda)))
        wrapper1 <- function(i) {
            tryCatch({
                fit = glmnet(t(Z[-i, ]), Z[i, ], family = link, 
                  lambda = lambda, standardize = standardize)
            }, error = function(e) {
                fit = glmnetEmpty(t(Z[-i, ]), lambda)
            })
            fit$beta = as.matrix(fit$beta)
            if (ncol(fit$beta) < length(lambda)) {
                tmp = matrix(0, nrow = nrow(fit$beta), ncol = length(lambda))
                tmp[, 1:ncol(fit$beta)] = fit$beta
                tmp[, ncol(fit$beta):length(lambda)] = fit$beta[, 
                  ncol(fit$beta)]
                fit$beta = tmp
            }
            if (i == 1) {
                ghat[i, 2:nrow(Z), ] = fit$beta
            }
            else if (i == nrow(Z)) {
                ghat[i, 1:(nrow(Z) - 1), ] = fit$beta
            }
            else {
                ghat[i, 1:(i - 1), ] = fit$beta[1:(i - 1), ]
                ghat[i, (i + 1):nrow(Z), ] = fit$beta[i:nrow(fit$beta), 
                  ]
            }
            return(ghat[i, , ])
        }
        if (parallel) {
            if (q == 0) {
                ghat2 = mclapply(1:nrow(Z), wrapper1)
                for (i in 1:nrow(Z)) {
                  ghat[i, , ] = ghat2[[i]]
                }
            }
            if (q != 0) {
                ghat2 = mclapply((p + 1):nrow(Z), wrapper1)
                for (i in (p + 1):nrow(Z)) {
                  ghat[i, , ] = ghat2[[i - p]]
                }
            }
            return(ghat)
        }
        if (parallel == F) {
            if (q == 0) {
                ghat2 = lapply(1:nrow(Z), wrapper1)
                for (i in 1:nrow(Z)) {
                  ghat[i, , ] = ghat2[[i]]
                }
            }
            if (q != 0) {
                ghat2 = lapply((p + 1):nrow(Z), wrapper1)
                for (i in (p + 1):nrow(Z)) {
                  ghat[i, , ] = ghat2[[i - p]]
                }
            }
            return(ghat)
        }
    }
    if (length(lambda) == 1) {
        ghat = matrix(0, nrow = nrow(Z), ncol = nrow(Z))
        if (parallel) {
            library(snowfall)
            sfInit(parallel = TRUE, cpus = nCpus)
            sfExport("X", local = T)
            sfExport("ghat", local = T)
            sfLibrary(glmnet)
            wrapper2 <- function(i) {
                tryCatch({
                  fit = glmnet(t(Z[-i, ]), Z[i, ], family = link, 
                    lambda = lambda, standardize = standardize)
                }, error = function(e) {
                  fit = glmnetEmpty(t(Z[-i, ]), lambda)
                })
                fit$beta = as.numeric(fit$beta)
                if (i == 1) {
                  ghat[i, 2:nrow(Z)] = fit$beta
                }
                else if (i == nrow(Z)) {
                  ghat[i, 1:(nrow(Z) - 1)] = fit$beta
                }
                else {
                  ghat[i, 1:(i - 1)] = fit$beta[1:(i - 1)]
                  ghat[i, (i + 1):nrow(Z)] = c(fit$beta[i:length(fit$beta)])
                }
                return(ghat[i, ])
            }
            if (q == 0) {
                sfExport("wrapper2")
                ghat = sfSapply(1:nrow(Z), wrapper2)
                sfStop()
            }
            if (q != 0) {
                sfExport("wrapper2")
                ghat = sfSapply((p + 1):nrow(Z), wrapper2)
                sfStop()
            }
            return(ghat)
        }
        if (parallel == F) {
            st = p + 1
            if (q == 0) {
                st = 1
            }
            for (i in st:nrow(Z)) {
                tryCatch({
                  fit = glmnet(t(Z[-i, ]), Z[i, ], family = link, 
                    lambda = lambda, standardize = standardize)
                }, error = function(e) {
                  fit = glmnetEmpty(t(Z[-i, ]), lambda)
                })
                fit$beta = as.numeric(fit$beta)
                if (i == 1) {
                  ghat[i, 2:nrow(Z)] = fit$beta
                }
                else if (i == nrow(Z)) {
                  ghat[i, 1:(nrow(Z) - 1)] = fit$beta
                }
                else {
                  ghat[i, 1:(i - 1)] = fit$beta[1:(i - 1)]
                  ghat[i, (i + 1):nrow(Z)] = c(fit$beta[i:length(fit$beta)])
                }
            }
            return(ghat)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
