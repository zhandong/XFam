glmpois <-
function(X, lambda, parallel=FALSE, nCpus = 4){
	
	if(length(lambda)>1){
		ghat = array(0,dim=c(nrow(X),nrow(X),length(lambda)))
		
		if(parallel){
			wrapper <- function(i){
				fit=glmnet(t(X[-i,]),X[i,],family="poisson",lambda= lambda)
				fit$beta=as.matrix(fit$beta)
				if(ncol(fit$beta)<length(lambda)){
					tmp = matrix(0,nrow = nrow(fit$beta),ncol = length(lambda))
					tmp[,1:ncol(fit$beta)]=fit$beta
					tmp[,ncol(fit$beta):length(lambda)] = fit$beta[,ncol(fit$beta)]
					fit$beta = tmp
				}
				
				if(i==1){
					ghat[i,2:nrow(X),]=fit$beta
				}else if(i==nrow(X)){
					ghat[i,1:(nrow(X)-1),]=fit$beta
				}else{
					ghat[i,1:(i-1),]=fit$beta[1:(i-1),]
					ghat[i,(i+1):nrow(X),]=fit$beta[i:nrow(fit$beta),]	
				}
				return(ghat[i,,])
			}

			library(multicore)
			ghat2=mclapply(1:nrow(X),wrapper)		
			
			for(i in 1:nrow(X)){
				ghat[i,,]=ghat2[[i]]
			}
			
			return(ghat)
		}
	
		if(parallel==F){
			wrapper <- function(i){
				#print(i)
				fit=glmnet(t(X[-i,]),X[i,],family="poisson",lambda= lambda)
				fit$beta=as.matrix(fit$beta)
				if(ncol(fit$beta)<length(lambda)){
					tmp = matrix(0,nrow = nrow(fit$beta),ncol = length(lambda))
					tmp[,1:ncol(fit$beta)]=fit$beta
					tmp[,ncol(fit$beta):length(lambda)] = fit$beta[,ncol(fit$beta)]
					fit$beta = tmp
				}
				
				if(i==1){
					ghat[i,2:nrow(X),]=fit$beta
				}else if(i==nrow(X)){
					ghat[i,1:(nrow(X)-1),]=fit$bet
				}else{
					ghat[i,1:(i-1),]=fit$beta[1:(i-1),]
					ghat[i,(i+1):nrow(X),]=fit$beta[i:nrow(fit$beta),]	
				}
				return(ghat[i,,])
			}
			ghat2=lapply(1:nrow(X),wrapper)	
			for(i in 1:nrow(X)){
				ghat[i,,]=ghat2[[i]]
			}
			return(ghat)
		}
	}
	
	if(length(lambda)==1){		
		ghat=matrix(0,nrow=nrow(X),ncol=nrow(X))
		
		if(parallel){
			library(snowfall)
			sfInit(cpus=nCpus)
			
			sfExport("X",local=T)
			sfExport("ghat",local=T)
			sfLibrary(glmnet)
		#-modify ghat
			wrapper <- function(i){
				fit=glmnet(t(X[-i,]),X[i,],family="poisson",lambda= lambda)
				fit$beta=as.numeric(fit$beta)
				if(i==1){
					ghat[i,2:nrow(X)]=fit$beta
				}else if(i==nrow(X)){
					ghat[i,1:(nrow(X)-1)]=fit$beta
				}else{
					ghat[i,1:(i-1)]=fit$beta[1:(i-1)]
					ghat[i,(i+1):nrow(X)]=c(fit$beta[i:length(fit$beta)])	
				}
				return(ghat[i,])
			}
			sfExport("wrapper")
			ghat=sfSapply(1:nrow(X),wrapper)	
			sfStop()
			return(ghat)
		}
		
# wooi question: should run this again if parallel==F?
		for(i in 1:nrow(X)){
			fit=glmnet(t(X[-i,]),X[i,],family="poisson",lambda= lambda)
			fit$beta=as.numeric(fit$beta)
			if(i==1){
				ghat[i,2:nrow(X)]=fit$beta
			}else if(i==nrow(X)){
				ghat[i,1:(nrow(X)-1)]=fit$beta
			}else{
				ghat[i,1:(i-1)]=fit$beta[1:(i-1)]
				ghat[i,(i+1):nrow(X)]=c(fit$beta[i:length(fit$beta)])	
			}
				
		}
		return(ghat)
	}
}
