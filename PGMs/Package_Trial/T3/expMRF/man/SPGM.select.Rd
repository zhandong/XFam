\name{SPGM.select}
\alias{SPGM.select}
\title{Log-Linear Graphical Model based on Pair-wise Sub-linear truncated Poisson Markov Network}
\description{
	Fitting the local Log-Linear Graphical Model based on pair-wise sublinear-truncated Poisson Markov Network. The network modeling algorithm is the same as LPGM.
}
\usage{
SPGM.select(X, R, R0 = 0, N = 100, beta = 0.05, lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{R}{lower-bound threshold value for the trunctation, has to be positive}
  \item{R0}{lower-bound threshold value for the trunctation, default to 0}
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

\seealso{
	\code{\link{GMS}},
	\code{\link{Bsublin}}, 
	\code{\link{LPGM.select}}
}
\examples{
library(PGM)
library(huge)
n = 200
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace=TRUE), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

range(simData)
spgm.path.all.p = SPGM.select(t(simData), 4, 2, nlams=20, N=10, beta=0.05, nCpus=2, parallel=TRUE)
spgm.path.all.p


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, R0 = 0, N = 100, beta = 0.05, lmin = 0.01, nlams = 20, 
    lambda.path = NULL, parallel = T, nCpus = 4) 
{
    require("huge")
    require("glmnet")
    if (R < 0) {
        cat("ERROR: Truncating threshold R should be positive. \n")
        ghat = NULL
        return(ghat)
    }
    Xorig <- X
    X <- round(Bsublin(X, R, R0))
    ghat <- LPGM.select(X, method = "SPGM", N = N, beta = beta, 
        lmin = lmin, nlams = nlams, lambda.path = lambda.path, 
        parallel = parallel, nCpus = nCpus)
    if (!is.null(ghat)) {
        ghat$call = match.call()
    }
    return(ghat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
