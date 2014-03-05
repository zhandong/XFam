\name{SPGM.select}
\alias{SPGM.select}
\title{Log-Linear Graphical Model based on Pair-wise Sub-linear truncated Poisson Markov Network}
\description{
	Fitting the local Log-Linear Graphical Model based on pair-wise sublinear-truncated Poisson Markov Network. The network modeling algorithm is the same as LPGM.
}
\usage{
SPGM.select(X, R, R0 = 0, N = 100, beta = 0.05, lmin = 0.01, nlams = 20, 
lambda.path = NULL, parallel = T, nCpus = 4)
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
	\code{\link{Bsublin}}, 
	\code{\link{LPGM.select}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, R, R0 = 0, N = 100, beta = 0.05, lmin = 0.01, nlams = 20, 
    lambda.path = NULL, parallel = T, nCpus = 4) 
{
    if (R < 0) {
        cat("ERROR: Truncating threshold R should be positive. \n")
        ghat = NULL
        return(ghat)
    }
    Xorig <- X
    X <- round(Bsublin(X, R, R0))
    return(LPGM.select(X, method = "SPGM", N = N, beta = beta, 
        lmin = lmin, nlams = nlams, lambda.path = lambda.path, 
        parallel = parallel, nCpus = nCpus))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
