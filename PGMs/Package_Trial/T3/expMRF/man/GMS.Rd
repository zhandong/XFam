\name{GMS}
\alias{GMS}
\title{
Local Log-linear Graphical Models
}
\description{
This class of objects is returned by various "GM" functions included in this expMRF package, to represent the fitted markov networks over the regularization paths.
Objects of this class have the print method to display the core information of the fitted models and plot method to plot the optimal markov network.
}
\arguments{
  \item{v}{vector of (nlams) variability measured from the stability selection}
  \item{lambda.path}{vector lambda used for regularization}
  \item{opt.lambda}{lambda value that gives the optimal network (network with maximum variability)}
  \item{network}{a list of pxp coefficient matrix along the regularization.}
  \item{opt.index}{index of the regularization value that gives the optimal network}
}
\details{
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

\seealso{
	\code{\link{LPGM.select}}, 
	\code{\link{SPGM.select}},
	\code{\link{TPGM.select}},
	\code{\link{LGGM.select}},
	\code{\link{LISM.select}}
}