LPGM.select <-
function(X,method="LPGM",N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=TRUE,nCpus=4){
	ghat <- LGM.select.generic(X, method=method, link="poisson", N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, nCpus=nCpus)
	if(!is.null(ghat)){
		ghat$call <- match.call()
	}
	
	return(ghat)
}
