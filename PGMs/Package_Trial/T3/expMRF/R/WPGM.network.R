WPGM.network <-
function(X,R,nlams,lmin=0.001,lambda=NULL, parallel=TRUE,ncores=4){
	if(is.null(lambda)){
		lmax = lambdaMax(t(X))
		lambda = exp(seq(log(lmax),log(lmin),l=nlams));	
	}
	if(nlams!= length(lambda)){
		print("nlams is not equal to lams")
	}
	ghat = c()
	if(nlams>0){
		ghat = array(0,dim=c(nrow(X),nrow(X),length(lambda)))
	}
	wrapper <- function(i){
		fit = WPGM.path.neighborhood(t(X[-i,]),X[i,],R,nlams,lambda=lambda,0)
		fit$beta=as.matrix(fit$Bmat)
		if(i==1){
			ghat[i,2:nrow(X),]=fit$beta
		}
		else if(i==nrow(X)){
			ghat[i,1:(nrow(X)-1),]=fit$beta
		}
		else{
			ghat[i,1:(i-1),]=fit$beta[1:(i-1),]
			ghat[i,(i+1):nrow(X),]=fit$beta[i:nrow(fit$beta),]	
		}
		return(ghat[i,,])
	}
	
	ghat2=c()
	if(parallel){
		library(multicore)
		ghat2=mclapply(1:nrow(X),wrapper,mc.cores=ncores)	
	}
	else{
		ghat2=lapply(1:nrow(X),wrapper)	
	}
	for(i in 1:nrow(X)){
		ghat[i,,]=ghat2[[i]]
	}

	ghat=lapply(1:nlams,function(r){return(ghat[,,r])})
	return(ghat)
}
