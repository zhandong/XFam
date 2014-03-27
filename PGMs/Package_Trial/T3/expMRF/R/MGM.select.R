MGM.select <-
function(X, Y, xlink="poisson", ylink="binomial", method="Both", N=100, beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL, parallel=TRUE,nCpus=4, standardize=TRUE){
	 
	if(is.null(Y)){
		#Z <- X
		#p <- nrow(Z)
		#q <- 0
		cat("ERROR: Y could not be empty. \n")
		ghat = NULL
		return(ghat)
	}
	
	if(!is.null(Y)){
		if(ncol(X) == ncol(Y)){
			Z <- rbind(X, Y)
			p = nrow(X)
			q = nrow(Y)
			if(is.null(rownames(X))){
				rownames(X) <- paste("X_V", 1:p, sep="")
			}
			if(is.null(rownames(Y))){
				rownames(Y) <- paste("Y_V", 1:q, sep="")
			}
		}
		else{
			cat("ERROR: But X and Y should have same number of observations. \n")
			ghat = NULL
			return(ghat)
		}
	}
	
	if(is.null(lambda.path) ){
		#lmax = myglmnet.max(Z)
		# transpose Z as our Z is pxn - theorectically max(X'X) where X is nxp
		tmp = t(t(Z))%*%t(Z)
		lmax = max(tmp[upper.tri(tmp)])
 		lambda.path = exp(seq(log(lmax),log(lmin),l=nlams))
	}
	#
	b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
	ghat=list()
	v=c()
	
	if(parallel == T){
		ghat.path=list()
		ghat.path$path=vector("list",length(lambda.path))
		
		for(i in 1:N){
			cat(paste("MGM: Conducting sampling ... in progress: ", floor(100*(i/N)), "%", collapse=""),"\r")
			flush.console()
			
			index = sample(1:ncol(Z),b,replace=F)
			
			#network.X <- glmGeneric(X[, index], Y[, index], link="binomial",lambda=lambda.path,parallel=parallel,nCpus=4, standardize=TRUE)
			#network.Y <- glmGeneric(Y[, index], X[, index], link="poisson",lambda=lambda.path,parallel=parallel,nCpus=4, standardize=TRUE)
			
			if (method == "Both"){
				network.X <- glmGeneric(X[, index], Y[, index], link=ylink, lambda=lambda.path, parallel=parallel, nCpus=4, standardize=TRUE)
				network.Y <- glmGeneric(Y[, index], X[, index], link=xlink, lambda=lambda.path, parallel=parallel, nCpus=4, standardize=TRUE)
			}
			
			if (method == "Right"){
				network.X <- glmGeneric(X[, index], Y[, index], link=ylink, lambda=lambda.path, parallel=parallel, nCpus=4, standardize=TRUE)
				network.Y <- glmGeneric(X[, index], NULL, link=xlink, lambda=lambda.path, parallel=parallel, nCpus=4, standardize=TRUE)
			}
			
			if (method == "Left"){
				network.X <- glmGeneric(Y[, index], NULL, link=ylink, lambda=lambda.path, parallel=parallel, nCpus=4, standardize=TRUE)
				network.Y <- glmGeneric(Y[, index], X[, index], link=xlink, lambda=lambda.path, parallel=parallel, nCpus=4, standardize=TRUE)
			}
			
			network.raw <- Merge.GraphXY(network.X, network.Y, rownames(X), rownames(Y), length(lambda.path), method)
			
			for(j in 1:length(lambda.path)){
				#tmp=ghat.path$raw[,,j]
				tmp = network.raw[[j]]
				tmp[abs(tmp)<1e-06]=0
				tmp[abs(tmp)>1e-06]=1
				diag(tmp)=0
				if(is.null(ghat.path$path[[j]])){
					ghat.path$path[[j]]=tmp;
				}else{
					ghat.path$path[[j]]=ghat.path$path[[j]]+tmp	
				}
			}
		}
		
		for(i in 1:length(lambda.path)){
			D=ghat.path$path[[i]]
			D=D/N
			D=2*D*(1-D)
			v=c(v,mean(D[upper.tri(D)]))	
		}
	}	

	#
	if(parallel == F){
		for( j in 1:length(lambda.path)){
			cat(paste("MGM: Conducting sampling ... in progress: ", floor(100*(j/length(lambda.path))), "%", collapse=""),"\r")
			flush.console()
			
			# Do regression on Y as the response - thus logistic regression
			D=matrix(0,nrow=nrow(Z),ncol=nrow(Z))
			for(i in 1:N){
				index = sample(1:ncol(Z),b,replace=F)
				#tmp=glmpois(X[,index],lambda.path[j],parallel=F)
				
				if (method == "Both"){
					#tmp1 = glmGeneric(X[,index], Y[,index], link="binomial", lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					tmp1 = glmGeneric(X[,index], Y[,index], link=ylink, lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					rownames(tmp1) <- colnames(tmp1) <- c(rownames(X), rownames(Y))
					
					#tmp2 = glmGeneric(Y[,index], X[,index], link="poisson", lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					tmp2 = glmGeneric(Y[,index], X[,index], link=xlink, lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					rownames(tmp2) <- colnames(tmp2) <- c(rownames(Y), rownames(X))
					
					tmp = tmp1
					tmp[rownames(X), colnames(tmp)] = tmp2[rownames(X), colnames(tmp)]
				}
				
				if (method == "Right"){
					tmp1 = glmGeneric(X[,index], Y[,index], link=ylink, lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					rownames(tmp1) <- colnames(tmp1) <- c(rownames(X), rownames(Y))
					
					tmp2 = glmGeneric(X[,index], NULL, link=xlink, lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					rownames(tmp2) <- colnames(tmp2) <- rownames(X)
					
					tmp = tmp1
					tmp[rownames(X), rownames(X)] = tmp2[rownames(X), rownames(X)]
				}
				
				if (method == "Left"){
					tmp1 = glmGeneric(Y[,index],NULL, link=ylink, lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					rownames(tmp1) <- colnames(tmp1) <- rownames(Y)
					
					tmp2 = glmGeneric(Y[,index], X[,index], link=xlink, lambda=lambda.path[j], parallel=parallel, nCpus=nCpus, standardize=TRUE)
					rownames(tmp2) <- colnames(tmp2) <- c(rownames(Y), rownames(X))
					
					tmp = tmp2
					tmp[rownames(Y), rownames(Y)] = tmp1[rownames(Y), rownames(Y)]
				}
				
				tmp[abs(tmp)<1e-06]=0
				tmp[abs(tmp)>1e-06]=1
				D=D+tmp
			}
			
			D=D/N
			D=2*D*(1-D)
			v=c(v,mean(D[upper.tri(D)]))			
		}
	}	
	#
	v=cummax(v)
	ghat$v=v
	ghat$lambda.path = lambda.path
	ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]	
	
	if (method == "Both"){
		#network.X <- glmGeneric(X,Y,link="binomial",lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		#network.Y <- glmGeneric(Y,X,link="poisson",lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		network.X <- glmGeneric(X,Y,link=ylink,lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		network.Y <- glmGeneric(Y,X,link=xlink,lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		#ghat$network <- Merge.GraphXY(network.X, network.Y, rownames(X), rownames(Y), length(lambda.path))
	}
	if (method == "Right"){
		#network.X <- glmGeneric(X,Y,link="binomial",lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		#network.Y <- glmGeneric(Y,X,link="poisson",lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		network.X <- glmGeneric(X,Y,link=ylink,lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		network.Y <- glmGeneric(X,NULL,link=xlink,lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		#ghat$network <- Merge.GraphXY(network.X, network.Y, rownames(X), rownames(Y), length(lambda.path))
	}
	if (method == "Left"){
		#network.X <- glmGeneric(X,Y,link="binomial",lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		#network.Y <- glmGeneric(Y,X,link="poisson",lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		network.X <- glmGeneric(Y,NULL,link=ylink,lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		network.Y <- glmGeneric(Y,X,link=xlink,lambda=lambda.path,parallel=parallel,nCpus=nCpus, standardize=TRUE)
		#ghat$network <- Merge.GraphXY(network.X, network.Y, rownames(X), rownames(Y), length(lambda.path))
	}
	
	ghat$network <- Merge.GraphXY(network.X, network.Y, rownames(X), rownames(Y), length(lambda.path), method)
	ghat$opt.index = which(v==max(v[v<beta]))
	ghat$call <- match.call()
	
	cat(paste("\n", "MGM Completed.", "\n", sep=""))
	
	class(ghat) <- "GMS"
	
	return(ghat)
}
