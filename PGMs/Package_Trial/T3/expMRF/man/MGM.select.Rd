\name{MGM.select}
\alias{MGM.select}
\title{Mixed-Graphical Model}
\description{
	Fitting the local Log-linear Graphical model based on pair-wise markov properties from a data matrix with combination of two data types (poisson and binary). 
}
\usage{
MGM.select(X, Y, xlink = "poisson", ylink = "binomial", method = "Both", N = 100, beta = 0.05, 
	lmin = 0.01, nlams = 20, lambda.path = NULL, parallel = TRUE, nCpus = 4, standardize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \ item{X}{a pxn data matrix.}
  \ item{Y}{a qxn data matrix.}
  \ item{xlink}{the distribution family for data X, default to "poisson"}
  \ item{ylink}{the distribution family for data Y, default to "binomial"}
  \ item{method}{specification of the variation of methods in inferring the network between two data types. Refer details section below for more discussion.}
  \ item{N}{number of iteration for stability selection, default to 100}
  \ item{beta}{threshold value on sparsity of the network to filter out dense network}
  \ item{lmin}{minimum lambda value, default to 0.01}
  \ item{nlams}{number of lambda for regularization}
  \ item{lambda.path}{vector lambda used for regularization}
  \ item{parallel}{logical value to indicate if the process should be run parallelly in multiple threads, default to TRUE}
  \ item{nCpus}{number of (maximum) cores to use for parallel execution, default to 4}
  \ item{standardize}{Logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is
standardize=TRUE. This parameter passed to glmnet for the parameter of the same name.}
}
\details{
To infer network of the data matrix Z (combined data from X and Y by row), three alternative methods are implemented: 
method = "Both", mixed MRF of Z (or X <-> Y); which are performing two CRF regressions: \code{X_i ~ X_{\ i} + Y} and \code{Y_i ~ X + Y_{\ i}}
method = "Right"; which are which are performing a MRF regression on X: \code{X_i ~ X_{\ i}} and CRF on Y \code{Y_i ~ X + Y_{\ i}}
method = "Left"; which are which are performing a CRF regression on X: \code{X_i ~ X_{\ i} + Y} and MRF on Y \code{Y_i ~ Y_{\ i}}
For each method, MRF or CRF are implemented separately and then only merge to form one final inferred network.

Also, the sample size of both data X and Y should be the same and they should be arranged into the same order.
}
\value{
an object of class GMS object will be returned, represents the modeled markov networks over the regularization path. See GMS for details.
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
	\code{\link{Merge.GraphXY}}
}
\examples{
n = 100
p = 100
q = 50

tData.X = matrix(rpois(n*p,1),n,p)
tData.Y = matrix(rbinom(n*q,1,0.5),n,q)

fixpath <- exp(seq(log(10),log(0.01),l=10))

mgm.nets <- MGM.select(t(tData.X),t(tData.Y),xlink="poisson", ylink="binomial", method="Both", N=10,beta=0.05, lmin = 0.01, nlams=length(fixpath), lambda.path=fixpath, parallel=F,standardize=TRUE)
mgm.nets.right <- MGM.select(t(tData.X),t(tData.Y),xlink="poisson", ylink="binomial", method="Right", N=10,beta=0.05, lmin = 0.01, nlams=length(fixpath), lambda.path=fixpath, parallel=F,standardize=TRUE)
mgm.nets.left <- MGM.select(t(tData.X),t(tData.Y),xlink="poisson", ylink="binomial", method="Left", N=10,beta=0.05, lmin = 0.01, nlams=length(fixpath), lambda.path=fixpath, parallel=F,standardize=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
