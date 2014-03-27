\name{LGM.select.generic}
\alias{LGM.select.generic}
\title{Generic local log-linear graphical model}
\description{
A generic function allow the fitting of the local log-Linear graphical model based on pair-wise markov properties using an efficient, parallel algorithm over a path of regularization parameters (lambda). 
This algorithm employs neighborhood selection in inferring network structure. Stability selection is used in selecting the optimal network.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LGM.select.generic(X, method = "LPGM", link = "poisson", N = 100, beta = 0.05, lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4)
}
\arguments{
  \item{X}{a pxn data matrix}
  \item{method}{specification of the variation of local log-linear graphical model to be fitted.  Default to pair-wise Poisson graphical models (LPGM). 
			Other methods allowed included "TPGM" fro truncated Poisson, "SPGM" for sublinear Poisson, "LGGM" for Gaussian, "LISM" for Ising model (binomial).}
  \item{link}{specification of the exponential family of the data matrix X.  Default to Poisson distribution ("poisson"). 
			Other links allowed "gaussian" for Gaussian, "binomial" for binary data.}
  \item{N}{number of iteration for stability selection, default to 100}
  \item{beta}{threshold value on sparsity of the network to filter out dense network, default to 0.05}
  \item{lmin}{minimum lambda value, default to 0.01}
  \item{nlams}{number of lambda for regularization, default to 20}
  \item{lambda.path}{vector lambda used for regularization, defaults to NULL}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to TRUE}
  \item{nCpus}{number of (maximum) cores to use for parallel execution, default to 4}
}
\details{
This is the generic function to implement the local log-Linear markov network (proposed in the first reference below).  
The core function of this funciton is to employ neighborhood selection to infer the network structure (\code{\link{glmGeneric}}) and STAR (stability selection) for various distribution families.  
}
\value{
an object of class GMS object will be returned, represents the modeled markov networks over the regularization path. See GMS for details.
}
\references{
G.I. Allen and Z. Liu, 2012, A Log-Linear Graphical Model for Inferring Genetic Networks from High-Throughput Sequencing Data, \emph{The IEEE International Conference on Bioinformatics and Biomedicine (BIBM 2012)}.

E. Yang, P.K. Ravikumar, G.I. Allen, and Z. Liu, 2012, Graphical Models via Generalized Linear Models, \emph{NIPS}, \bold{vol. 25}, pp. 1367--1375.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{GMS}},
	\code{\link{glmGeneric}}, 
	\code{\link{myglmnet.max}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, method = "LPGM", link = "poisson", N = 100, beta = 0.05, 
    lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = T, 
    nCpus = 4) 
{
    require("huge")
    require("glmnet")
    if (is.null(lambda.path)) {
        lmax = myglmnet.max(X, link = link)
        lambda.path = exp(seq(log(lmax), log(lmin), l = nlams))
    }
    if (parallel == T) {
        b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
        ghat = list()
        ghat.path = list()
        ghat.path$path = vector("list", length(lambda.path))
        v = c()
        for (i in 1:N) {
            cat(paste(method, ": Conducting sampling ... in progress: ", 
                floor(100 * (i/N)), "\%", collapse = ""), "\r")
            flush.console()
            glmpois.good <- 1
            while (glmpois.good) {
                good <- 1
                while (good) {
                  index = sample(1:ncol(X), b, replace = F)
                  if (sum(apply(X[, index], 1, function(x) length(unique(x)) == 
                    1)) == 0) {
                    good <- 0
                  }
                }
                tryCatch({
                  ghat.path$raw = glmGeneric(X[, index], NULL, 
                    link = link, lambda = lambda.path, parallel = T, 
                    nCpus = nCpus)
                  glmpois.good <- 0
                }, error = function(e) {
                  cat("glmnet returns empty model. Try again.")
                })
            }
            for (j in 1:length(lambda.path)) {
                tmp = ghat.path$raw[, , j]
                tmp[abs(tmp) < 1e-06] = 0
                tmp[abs(tmp) > 1e-06] = 1
                diag(tmp) = 0
                if (is.null(ghat.path$path[[j]])) {
                  ghat.path$path[[j]] = tmp
                }
                else {
                  ghat.path$path[[j]] = ghat.path$path[[j]] + 
                    tmp
                }
            }
        }
        for (i in 1:length(lambda.path)) {
            D = ghat.path$path[[i]]
            D = D/N
            D = 2 * D * (1 - D)
            v = c(v, mean(D[upper.tri(D)]))
        }
        v = cummax(v)
        ghat$v = v
        ghat$lambda.path = lambda.path
        ghat$opt.lambda = lambda.path[which(v == max(v[v < beta]))]
        ghat$network = glmGeneric(X, NULL, link = link, lambda = lambda.path, 
            parallel = T, nCpus = nCpus)
        ghat$network = lapply(1:nlams, function(r) {
            return(ghat$network[, , r])
        })
        ghat$opt.index = which(v == max(v[v < beta]))
        ghat$call <- match.call()
        cat(paste("\n", method, " Completed.", "\n", sep = ""))
        class(ghat) <- "GMS"
        return(ghat)
    }
    if (parallel == F) {
        b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
        ghat = list()
        v = c()
        for (j in 1:length(lambda.path)) {
            cat(paste(method, ": Conducting sampling ... in progress: ", 
                floor(100 * (j/length(lambda.path))), "\%", collapse = ""), 
                "\r")
            flush.console()
            D = matrix(0, nrow = nrow(X), ncol = nrow(X))
            for (i in 1:N) {
                glmpois.good <- 1
                while (glmpois.good) {
                  good <- 1
                  while (good) {
                    index = sample(1:ncol(X), b, replace = F)
                    if (sum(apply(X[, index], 1, function(x) length(unique(x)) == 
                      1)) == 0) {
                      good <- 0
                    }
                  }
                  tryCatch({
                    tmp = glmGeneric(X[, index], NULL, link = link, 
                      lambda = lambda.path[j], parallel = F)
                    glmpois.good <- 0
                  }, error = function(e) {
                    cat("glmnet returns empty model. Try again.\n")
                  })
                }
                tmp[abs(tmp) < 1e-06] = 0
                tmp[abs(tmp) > 1e-06] = 1
                D = D + tmp
            }
            D = D/N
            D = 2 * D * (1 - D)
            v = c(v, mean(D[upper.tri(D)]))
        }
        v = cummax(v)
        ghat$v = v
        ghat$lambda.path = lambda.path
        ghat$opt.lambda = lambda.path[which(v == max(v[v < beta]))]
        ghat$network = glmGeneric(X, NULL, link = link, lambda = lambda.path, 
            parallel = parallel, nCpus = nCpus)
        ghat$network = lapply(1:nlams, function(r) {
            return(ghat$network[, , r])
        })
        ghat$opt.index = which(v == max(v[v < beta]))
        ghat$call <- match.call()
        cat(paste("\n", method, " Completed.", "\n", sep = ""))
        class(ghat) <- "GMS"
        return(ghat)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
