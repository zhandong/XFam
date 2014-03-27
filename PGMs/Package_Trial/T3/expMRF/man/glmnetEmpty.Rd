\name{glmnetEmpty}
\alias{glmnetEmpty}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Empty local log-linear graphical model.
}
\description{
This function will return an null (empty) local log-linear graphical model.
}
\usage{
glmnetEmpty(X, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{lambda}{lambda vector used for regularization}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{fit}{a matrix of coefficients with zero values, in the dimension of (number of lambda) x p (number of variable) }
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, lambda) 
{
    fit = list()
    fit$a0 <- rep(0, length(lambda))
    fit$lambda <- lambda
    fit$df <- 0
    fit$dim <- c(ncol(X), length(lambda))
    fit$beta <- Matrix(0, ncol(X), length(lambda))
    rownames(fit$beta) <- colnames(X)
    colnames(fit$beta) <- paste("s", 0:(length(lambda) - 1), 
        sep = "")
    return(fit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
