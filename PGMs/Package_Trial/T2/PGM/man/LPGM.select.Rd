\name{LPGM.select}
\alias{LPGM.select}
\title{Log-Linear Graphical Model based on Pair-wise Poisson Markov Network}
\description{
	Fitting the local Log-Linear Graphical Model based on pair-wise Poisson Markov Network using an efficient, parallel algorithm named Poisson Graphical Lasso over a path of regularization parameters (lambda). This algorithm employs neighborhood selection in inferring network structure. Stability selection is used in selecting the optimal network.
}
\usage{
LPGM.select(X, method = "LPGM", N = 100, beta = 0.05, lmin = 0.01, nlams = 20, 
lambda.path = NULL, parallel = T, nCpus = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{method}{specification of the variation of log-linear Poisson-based graphical model (LPGM), default to "LPGM". 
  Other two methods allowed are truncated poisson graphical model (TPGM) and sub-linear poisson graphical model (SPGM).
}
  \item{N}{number of iteration for stability selection, default to 100}
  \item{beta}{threshold value on sparsity of the network to filter out dense network}
  \item{lmin}{minimum lambda value, default to 0.01}
  \item{nlams}{number of lambda for regularization}
  \item{lambda.path}{vector lambda used for regularization}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to TRUE}
  \item{nCpus}{number of (maximum) cores to use for parallel execution, default to 4}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of five elements:
  \item{v}{vector of (nlams) variability measured from the stability selection}
  \item{lambda.path}{vector lambda used for regularization}
  \item{opt.lambda}{lambda value that gives the optimal network (network with maximum variability)}
  \item{network}{a list of pxp coefficient matrix along the regularization.}
  \item{opt.index}{index of the regularization value that gives the optimal network}
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
	\code{\link{myglmnet.max}}, 
	\code{\link{glmpois}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(PGM)
library(huge)
n = 200
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace=TRUE), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

#-# Run LPGM
lpgm.path.all = LPGM.select(t(simData), nlams=20, N=10, beta=0.05, parallel=FALSE)
str(lpgm.path.all)

## The function is currently defined as
function (X, method = "LPGM", N = 100, beta = 0.05, lmin = 0.01, 
    nlams = 20, lambda.path = NULL, parallel = T, nCpus = 4) 
{
    if (is.null(lambda.path)) {
        lmax = myglmnet.max(X)
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
                  ghat.path$raw = glmpois(X[, index], lambda = lambda.path, 
                    parallel = T, nCpus = nCpus)
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
        ghat$network = glmpois(X, lambda = lambda.path, parallel = T, 
            nCpus = nCpus)
        ghat$network = lapply(1:nlams, function(r) {
            return(ghat$network[, , r])
        })
        ghat$opt.index = which(v == max(v[v < beta]))
        cat(paste("\n", method, " Completed.", "\n", sep = ""))
        return(ghat)
    }
    if (parallel == F) {
        b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
        ghat = list()
        v = c()
        for (j in 1:length(lambda.path)) {
            cat(paste(method, ": Conducting sampling ... in progress: ", 
                floor(100 * (i/N)), "\%", collapse = ""), "\r")
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
                    tmp = glmpois(X[, index], lambda = lambda.path[j], 
                      parallel = F)
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
        ghat$network = glmpois(X, lambda = lambda.path, parallel = parallel, 
            nCpus = nCpus)
        ghat$network = lapply(1:nlams, function(r) {
            return(ghat$network[, , r])
        })
        ghat$opt.index = which(v == max(v[v < beta]))
        cat(paste("\n", method, " Completed.", "\n", sep = ""))
        return(ghat)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
