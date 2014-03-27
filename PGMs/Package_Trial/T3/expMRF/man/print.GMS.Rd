\name{print.GMS}
\alias{print.GMS}
\title{
	Print the GMS object.
}
\description{
	Default function to print the GMS object to the standard output.
}
\usage{
print.GMS(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a GMS object}
}
\details{
This is the default print function to display information of a GMS object to the standard output in readable format.
It will print the information of the function called to create the object, the index of the optimal network, the whole regularization path, estimated variabilities over the regularization path, 
and list of inferred markov networks over the whole regularization path.
Refer to \code{\link{GMS}} for details on GMS object.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
	\code{\link{GMS}}
}
\examples{
n = 100
p = 50
gdata = huge.generator(n,d=p, graph="scale-free",v=0.1,u=0.01)
smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace =T), nrow = nrow(gdata$theta) )
simData = WPGMSim(n,p,R=10, alpha = rep(0,p), Theta = 0.1*as.matrix(gdata$theta)*smatrix, maxit = 100 )

# Run LPGM
lpgm.path.all.p = LPGM.select(t(simData), nlams=10, N=10, beta=0.05, nCpus=2, parallel=T)

# Call this default print function
lpgm.path.all.p

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
