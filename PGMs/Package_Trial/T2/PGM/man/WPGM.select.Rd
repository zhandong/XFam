\name{WPGM.select}
\alias{WPGM.select}
\title{
	Winsorized Poisson Graphical Model (WPGM)
}
\description{
	Fitting the WPGM using efficient, parallel algorithm named Poisson Graphical Lasso.
	This algorithm employs neighborhood selection to infer network structure. 
	Stability selection method "star" was used in selecting the optimal network.
}
\usage{
WPGM.select(X, R=max(X), N=100, beta=0.05, lmin=0.0001, nlams=20, 
lambda.path=NULL, parallel=F, ncores = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{pxn data matrix}
  \item{R}{threshold value for truncating, default to be the maximum of value of the input data matrix}
  \item{N}{number of iteration for stability selection, default to 100}
  \item{beta}{threshold value on sparsity of the network to filter out dense network}
  \item{lmin}{minimum lambda value, default to 0.0001}
  \item{nlams}{number of lambda for regularization}
  \item{lambda.path}{vector lambda used for regularization}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to FALSE}
  \item{ncores}{number of (maximum) cores to use for parallel execution, default to 4}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R = max(X), method = "star", N = 100, beta = 0.05, 
    lambda.path = NULL, nlams = 20, ncores = 4, parallel = F) 
{
    if (is.null(lambda.path)) {
        lmax = lambdaMax(t(X))
        lambda.path = exp(seq(log(lmax), log(1e-04), l = nlams))
    }
    b = min(c(10 * sqrt(ncol(X)), 0.8 * ncol(X)))
    ghat = list()
    ghat.path = list()
    ghat.path$path = vector("list", length(lambda.path))
    v = c()
    for (i in 1:N) {
        cat(paste("WPGM: Conducting sampling ... in progress: ", 
            floor(100 * (i/N)), "\%", collapse = ""), "\r")
        flush.console()
        index = sample(1:ncol(X), b, replace = F)
        ghat.path$raw = WPGM.network(X[, index], R, nlams = length(lambda.path), 
            lambda = lambda.path, parallel = parallel, ncores = ncores)
        for (j in 1:length(lambda.path)) {
            tmp = ghat.path$raw[[j]]
            tmp[abs(tmp) < 1e-06] = 0
            tmp[abs(tmp) > 1e-06] = 1
            diag(tmp) = 0
            if (is.null(ghat.path$path[[j]])) {
                ghat.path$path[[j]] = tmp
            }
            else {
                ghat.path$path[[j]] = ghat.path$path[[j]] + tmp
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
    ghat$network = WPGM.network(X, R, nlams = length(lambda.path), 
        lambda = lambda.path, parallel = T)
    ghat$opt.index = which(v == max(v[v < beta]))
    cat("\nWPGM Completed. \n")
    return(ghat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
