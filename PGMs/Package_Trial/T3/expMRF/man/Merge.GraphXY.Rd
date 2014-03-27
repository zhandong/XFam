\name{Merge.GraphXY}
\alias{Merge.GraphXY}
\title{ Merge Graphs}
\description{
	Internal function to merge the two inferred markov networks produced by mixed-grpahical model (MGM.select function)
}
\usage{
Merge.GraphXY(network.X, network.Y, namesX, namesY, nlams, method = "Both")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network.X}{the first network object inferred}
  \item{network.Y}{the second network object inferred}
  \item{namesX}{the (p) variable names from data matrix X}
  \item{namesY}{the (q) variable naems from data matrix Y}
  \item{nlams}{number of lambda for regularization}
  \item{method}{specification of the variation of methods in inferring the network, default to "Both", 
  two other methods allowed are "Right" and "Left". Refer to \code{\link{MGM.select}} for details.
%%     ~~Describe \code{method} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
	A list of nlams merged networks with dimension (p+q) x (p+q); given X has p variables and Y has q variables.
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
	\code{\link{MGM.select}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (network.X, network.Y, namesX, namesY, nlams, method = "Both") 
{
    if (method == "Both") {
        results <- lapply(1:nlams, function(r) {
            tmp1 <- network.X[, , r]
            rownames(tmp1) <- colnames(tmp1) <- c(namesX, namesY)
            tmp2 <- network.Y[, , r]
            rownames(tmp2) <- colnames(tmp2) <- c(namesY, namesX)
            tmp = tmp1
            tmp[namesX, colnames(tmp)] = tmp2[namesX, colnames(tmp)]
            return(tmp)
        })
    }
    if (method == "Right") {
        results <- lapply(1:nlams, function(r) {
            tmp1 <- network.X[, , r]
            rownames(tmp1) <- colnames(tmp1) <- c(namesX, namesY)
            tmp2 <- network.Y[, , r]
            rownames(tmp2) <- colnames(tmp2) <- namesX
            tmp = tmp1
            tmp[namesX, namesX] = tmp2[namesX, namesX]
            return(tmp)
        })
    }
    if (method == "Left") {
        results <- lapply(1:nlams, function(r) {
            tmp1 <- network.X[, , r]
            rownames(tmp1) <- colnames(tmp1) <- namesY
            tmp2 <- network.Y[, , r]
            rownames(tmp2) <- colnames(tmp2) <- c(namesY, namesX)
            tmp = tmp2
            tmp[namesY, namesY] = tmp1[namesY, namesY]
            return(tmp)
        })
    }
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
