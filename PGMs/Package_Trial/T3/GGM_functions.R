#################
# Gaussian graphical model - local neighborhood seleciton
# depends on huge and glmnet packages
require('huge')
require('glmnet')

LGGM.select <- function(X,method="LGGM",N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=T,nCpus=4){
	ghat <- LGM.select.generic(X, method=method, link="gaussian", N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, nCpus=nCpus)
	if(!is.null(ghat)){
		ghat$call <- match.call()
	}
	
	return(ghat)
}
