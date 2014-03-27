\name{expMRF-package}
\alias{expMRF-package}
\alias{expMRF}
\docType{package}
\title{
Markov random field (MRF, or Markov network) for exponential family.
}
\description{
This package includes multiple functions to implement the local Log-Linear Graphical Model based on pair-wise markov network using efficient, parallel algorithm.
	This algorithm employs neighborhood selection to infer network structure. Stability selection is used in selecting the optimal network.
}
\details{
\tabular{ll}{
Package: \tab expMRF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-25\cr
License: \tab What license is it under?\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Ying-Wooi Wan, Zhandong Liu
Maintainer: Ying-Wooi Wan <yingwoow@bcm.edu>
~~ The author and/or maintainer of the package ~~
}
\references{
G.I. Allen and Z. Liu, 2012, A Log-Linear Graphical Model for Inferring Genetic Networks from High-Throughput Sequencing Data, \emph{The IEEE International Conference on Bioinformatics and Biomedicine (BIBM 2012)}.

E. Yang, P.K. Ravikumar, G.I. Allen, and Z. Liu, 2012, Graphical Models via Generalized Linear Models, \emph{NIPS}, \bold{vol. 25}, pp. 1367--1375.

}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
library(expMRF)

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
