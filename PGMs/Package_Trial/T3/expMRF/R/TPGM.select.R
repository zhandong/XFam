TPGM.select <-
function(X, R, N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=TRUE,nCpus=4){
	require('huge')
	require('glmnet')
	
	if (R < 0){
		cat("ERROR: Truncating threshold R should be positive. \n")
		ghat = NULL
		return(ghat)
	}
	
	# Should we check for other condition on R??
	
	# Transform the matrix with values truncated at R
	Xorig <- X
	X[X > R] <- R
	
	ghat <- LPGM.select(X, method="TPGM", N=N, beta=beta, lmin=lmin, nlams=nlams, lambda.path=lambda.path, parallel=parallel, nCpus=nCpus)
	if(!is.null(ghat)){
		ghat$call = match.call()
	}
	
	return(ghat)
}
