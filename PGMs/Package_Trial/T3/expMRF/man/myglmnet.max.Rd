\name{myglmnet.max}
\alias{myglmnet.max}
\title{Maximum lambda from binary search}
\description{
Obtain the regularization paramter lambda through binary search between zero to the maximum of X'X, in search for the smallest value that gives a null graphical model (empty network).
}
\usage{
myglmnet.max(X, link = "poisson", delta = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{link}{specification of the exponential family of the data matrix X.  Default to Poisson distribution ("poisson"). 
			Other links allowed "gaussian" for Gaussian, "binomial" for binary data.}
  \item{delta}{shift-size for the binary search, default to 0.01}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
numeric value for regularization parameter that will return a null model: the maximum lambda.
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
	\code{\link{lambdaMax}},
	\code{\link{glmGeneric}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(PGM)
library(huge)
library(glmnet)
n = 200
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace =TRUE), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

lmax = myglmnet.max(t(simData))


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, link = "poisson", delta = 0.01) 
{
    minlambda = 0
    maxlambda = lambdaMax(t(X))
    while (1) {
        mid = (minlambda + maxlambda)/2
        tmp = glmGeneric(X, NULL, link = link, lambda = mid)
        tmp[abs(tmp) < 1e-06] = 0
        tmp[abs(tmp) > 1e-06] = 1
        if (sum(tmp) > 0) {
            minlambda = mid + delta
        }
        else {
            maxlambda = mid - delta
        }
        if (abs(maxlambda - minlambda) < delta) {
            return(mid)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
