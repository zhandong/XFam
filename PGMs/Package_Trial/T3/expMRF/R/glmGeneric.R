glmGeneric <-
function(X, Y=NULL, link, lambda, parallel=FALSE, nCpus = 4, standardize=TRUE){
	if(is.null(Y)){
		Z <- X
		p <- nrow(Z)
		q <- 0
	}
	
	if(!is.null(Y)){
		if(ncol(X) == ncol(Y)){
			Z <- rbind(X, Y)
			p = nrow(X)
			q = nrow(Y)
		}
	}
	

	if(length(lambda)>1){
		#ghat = array(0,dim=c(nrow(X),nrow(X),length(lambda)))
		ghat = array(0,dim=c(nrow(Z),nrow(Z),length(lambda)))

		wrapper1 <- function(i){
			tryCatch(
				{
					fit = glmnet(t(Z[-i,]),Z[i,],family=link,lambda=lambda,standardize=standardize)
				},
				error = function(e) {
					fit = glmnetEmpty(t(Z[-i,]), lambda)
				})
			
			fit$beta=as.matrix(fit$beta)
			
			if(ncol(fit$beta)<length(lambda)){
				tmp = matrix(0,nrow = nrow(fit$beta),ncol = length(lambda))
				tmp[,1:ncol(fit$beta)]=fit$beta
				tmp[,ncol(fit$beta):length(lambda)] = fit$beta[,ncol(fit$beta)]
				fit$beta = tmp
			}
			
			if(i==1){
				ghat[i,2:nrow(Z),]=fit$beta
			}else if(i==nrow(Z)){
				ghat[i,1:(nrow(Z)-1),]=fit$beta
			}else{
				ghat[i,1:(i-1),]=fit$beta[1:(i-1),]
				ghat[i,(i+1):nrow(Z),]=fit$beta[i:nrow(fit$beta),]	
			}
			return(ghat[i,,])
		}

		if(parallel){
			if(q == 0){
				ghat2 = mclapply(1:nrow(Z),wrapper1)	
				for(i in 1:nrow(Z)){
					ghat[i,,]=ghat2[[i]]
				}
			}
			if(q != 0){
				ghat2 = mclapply((p+1):nrow(Z),wrapper1)
	# wooi: this might need to change ...because 1:p is nothing?!
				for(i in (p+1):nrow(Z)){
					ghat[i,,]=ghat2[[i-p]]
				}
			}
			return(ghat)
		}
	
		if(parallel==F)
		{
			if(q == 0){
				ghat2=lapply(1:nrow(Z),wrapper1)	
				for(i in 1:nrow(Z)){
					ghat[i,,]=ghat2[[i]]
				}
			}
			if(q != 0){
				ghat2=lapply((p+1):nrow(Z),wrapper1)
	# wooi: this might need to change ...because 1:p is nothing?!
				for(i in  (p+1):nrow(Z)){
					ghat[i,,]=ghat2[[i-p]]
				}
			}
			return(ghat)
		}
	}
	
	if(length(lambda) ==1){		
		ghat=matrix(0,nrow=nrow(Z),ncol=nrow(Z))
		
		if(parallel){
			library(snowfall)
			sfInit(parallel=TRUE, cpus=nCpus)
			
			sfExport("X",local=T)
			sfExport("ghat",local=T)
			sfLibrary(glmnet)
			
			wrapper2 <- function(i){
				tryCatch(
					{
						fit = glmnet(t(Z[-i,]),Z[i,],family=link,lambda= lambda,standardize=standardize)
					},
					error = function(e) {
						fit = glmnetEmpty(t(Z[-i,]), lambda)
					}
				)

				fit$beta=as.numeric(fit$beta)
				if(i==1){
					ghat[i,2:nrow(Z)]=fit$beta
				}else if(i==nrow(Z)){
					ghat[i,1:(nrow(Z)-1)]=fit$beta
				}else{
					ghat[i,1:(i-1)]=fit$beta[1:(i-1)]
					ghat[i,(i+1):nrow(Z)]=c(fit$beta[i:length(fit$beta)])	
				}
				return(ghat[i,])
			}
			
			if(q == 0){
				sfExport("wrapper2")
				ghat=sfSapply(1:nrow(Z),wrapper2)	
				sfStop()				
			}
			if(q != 0){
				sfExport("wrapper2")
				ghat=sfSapply((p+1):nrow(Z),wrapper2)	
				sfStop()	
			}
			return(ghat)
		}
		# wooi question: should run this again if parallel==T?
		if(parallel == F){
			st = p+1
			if(q == 0){
				st = 1
			}
			for(i in st:nrow(Z)){
				#cat(i)
				tryCatch(
					{
						fit = glmnet(t(Z[-i,]),Z[i,],family=link,lambda= lambda,standardize=standardize)
					},
					error = function(e) {
						fit = glmnetEmpty(t(Z[-i,]), lambda)
					}
				)
				fit$beta=as.numeric(fit$beta)
				if(i==1){
					ghat[i,2:nrow(Z)]=fit$beta
				}else if(i==nrow(Z)){
					ghat[i,1:(nrow(Z)-1)]=fit$beta
				}else{
					ghat[i,1:(i-1)]=fit$beta[1:(i-1)]
					ghat[i,(i+1):nrow(Z)]=c(fit$beta[i:length(fit$beta)])	
				}
			}
			return(ghat)
		}
	}
}
