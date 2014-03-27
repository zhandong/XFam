\name{plot.GMS}
\alias{plot.GMS}
\title{
	Plot GMS object.
}
\description{
	Default function to plot the optimal network of the GMS object (optimal markov network over the regularization path)
}
\usage{
plot.GMS(x, fn = "", th = 1e-06, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a GMS object}
  \item{fn}{file name to save the network plot. Default to be an empty string, 
  so the network is plotted to the standard output (screen). NOTE: if a filename is specificed, it should be filename for PDF file.}
  \item{th}{numeric value, default to 1e-06. To specify the threshold if the estimated coefficient between two variables is to be considered connected.}
  %\item{\dots}{other generic properties from the GMS object.}
}
\details{
	This is the default plotting function for GMS objects (markov networks inferred over a regularization path). Refer to \code{\link{GMS}} for details on GMS object.
	The function will plot the optimal network on the screen by default.  However, given a filename, the plot will be saved to a PDF file.
	The optimal netowrk will be plotted in force-directed layout (layout.kamada.kawai with default parameters implemented in igraph package).
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
lpgm.path.all.p

plot(lpgm.path.all.p, fn="lpgm.opt.net.pdf")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
