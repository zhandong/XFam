\name{LGGM.select}
\alias{LGGM.select}
\title{Local log-linear graphical model based on pair-wise Gaussian markov network.}
\description{
Fitting the local log-Linear graphical model based on pair-wise Gaussian markov network using an efficient, parallel algorithm over a path of regularization parameters (lambda). 
This algorithm employs neighborhood selection in inferring network structure. Stability selection is used in selecting the optimal network.
}
\usage{
LGGM.select(X, method = "LGGM", N = 100, beta = 0.05, lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a pxn data matrix}
  \item{method}{specification of the variation of local log-linear Gaussian-based graphical model (LGGM), default to "LGGM".}
  \item{N}{number of iteration for stability selection, default to 100}
  \item{beta}{threshold value on sparsity of the network to filter out dense network}
  \item{lmin}{minimum lambda value, default to 0.01}
  \item{nlams}{number of lambda for regularization}
  \item{lambda.path}{vector lambda used for regularization}
  \item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to TRUE}
  \item{nCpus}{number of (maximum) cores to use for parallel execution, default to 4}
}
\details{
This function is more of the interface to model the local log-linear Gaussian markov network.  Refer to \code{\link{LGM.select.generic}} for details in the model fitting.
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
	\code{\link{LGM.select.generic}}, 
	\code{\link{GMS}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	library(expMRF)
	n = 100
	p = 50
	tData.G = matrix(rnorm(n*p),n,p)
	lggm.path.all.p = LGGM.select(t(tData.G), nlams=10, N=10, beta=0.05, nCpus=2, parallel=TRUE)
	lggm.path.all.p
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
