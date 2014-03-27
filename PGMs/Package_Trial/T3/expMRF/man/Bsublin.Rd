\name{Bsublin}
\alias{Bsublin}
\title{Sublinear function}
\description{
	Transform the value of a data matrix (X) by a sub-linear function
}
\usage{
Bsublin(X, R, R0 = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a data matrix}
  \item{R}{upper-bound threshold value. Note: R should be great than 0}
  \item{R0}{lower-bound threshold value, default to 0}
}
\details{
Given two threshold values R and R0, s.t R > 0, R0 > 0 and R > R0

Each element \code{x} in X is transformed as follows:
	\code{\tabular{ll}{
		x = x,									\tab if x <= R0 \cr
		x = (-x^{2}+ 2*R*x - R0^{2})/(2 * (R - R0)), 	\tab if R0 < x <= R \cr
		x = (R + R0)/2,								\tab if x > R  \cr
	}}
}
\value{
	\item{Bx}{the transformed data matrix, of the same dimension as original data matrix X}
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
function (X, R, R0 = 0) 
{
    Bx = X
    Bx[X > R] = (R + R0)/2
    ind = X > R0 & X <= R
    Bx[ind] = (-X[ind]^2 + 2 * R * X[ind] - R0^2)/(2 * (R - R0))
    return(Bx)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
