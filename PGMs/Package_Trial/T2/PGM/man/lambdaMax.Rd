\name{lambdaMax}
\alias{lambdaMax}
\title{Maximum lambda}
\description{Compute the maximum lambda}
\usage{
lambdaMax(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{nxp data matrix}
}
\details{Largest value for regularization (maximum lambda), which is the maximum element from X'X}
\value{an integer value}
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
function (X) 
{
    tmp = t(X) \%*\% X
    return(max(tmp[upper.tri(tmp)]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
