\name{LISM.select}
\alias{LISM.select}
\title{Local log-linear graphical model based on Ising model.}
\description{
Fitting the Ising model through the efficient and parallel local log-linear graphical model algorithm over a path of regularization parameters. Stability selection is used in selecting the optimal network.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LISM.select(X, method = "LISM", N = 100, beta = 0.05, lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{method}{specification of the variation of local Ising-model (LISM), default to "LISM".}
  \item{N}{number of iteration for stability selection, default to 100}
  \item{beta}{threshold value on sparsity of the network to filter out dense network}
  \item{lmin}{minimum lambda value, default to 0.01}
  \item{nlams}{number of lambda for regularization}
  \item{lambda.path}{vector lambda used for regularization}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to TRUE}
  \item{nCpus}{number of (maximum) cores to use for parallel execution, default to 4}
}
\details{
This function is more of the interface to model the Ising model.  Refer to \code{\link{LGM.select.generic}} for details in the model fitting.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an object of class GMS object will be returned, represents the modeled markov networks over the regularization path. See GMS for details.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
	\code{\link{LGM.select.generic}}, 
	\code{\link{GMS}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(expMRF)

n = 100
p = 50
tData.Y = matrix(rbinom(n*p,1,0.5),n,p)
lism.path.all.p = LISM.select(t(tData.Y), nlams=10, N=10, beta=0.05, nCpus=2, parallel=TRUE)
lism.path.all.p


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, method = "LISM", N = 100, beta = 0.05, lmin = 0.01, 
    nlams = 20, lambda.path = NULL, parallel = T, nCpus = 4) 
{
    ghat <- LGM.select.generic(X, method = method, link = "binomial", 
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
