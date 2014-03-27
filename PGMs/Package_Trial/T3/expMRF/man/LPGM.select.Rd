\name{LPGM.select}
\alias{LPGM.select}
\title{Log-Linear Graphical Model based on Pair-wise Poisson Markov Network}
\description{
Fitting the local Log-Linear Graphical Model based on pair-wise Poisson Markov Network using an efficient, parallel algorithm named Poisson Graphical Lasso over a path of regularization parameters (lambda). This algorithm employs neighborhood selection in inferring network structure. Stability selection is used in selecting the optimal network.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LPGM.select(X, method = "LPGM", N = 100, beta = 0.05, lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4)
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
This function is more of the interface to model the local log-linear Poisson markov network.  Refer to \code{\link{LGM.select.generic}} for details in the model fitting.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{LGM.select.generic}}, 
	\code{\link{GMS}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(huge)
library(expMRF)
n = 200
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace=TRUE), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

#-# Run LPGM
lpgm.path.all = LPGM.select(t(simData), nlams=20, N=10, beta=0.05, parallel=FALSE)
str(lpgm.path.all)

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, method = "LPGM", N = 100, beta = 0.05, lmin = 0.01, 
    nlams = 20, lambda.path = NULL, parallel = T, nCpus = 4) 
{
    ghat <- LGM.select.generic(X, method = method, link = "poisson", 
        N = N, beta = beta, lmin = lmin, nlams = nlams, lambda.path = lambda.path, 
        parallel = parallel, nCpus = nCpus)
    if (!is.null(ghat)) {
        ghat$call <- match.call()
    }
    return(ghat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
