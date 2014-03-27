library(multicore)
require('huge')
require('glmnet')

# Reference : http://igraph.sourceforge.net/documentation.html


GMS <- function(x, ...) UseMethod("GMS")

# Default function to print the 'ghat' object 
print.GMS <- function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	
	cat(paste("\nOptimal lambda index: ", x$opt.index, "\n", sep=""))
	
	#cat(paste("\nLamba path (length: ", length(x$lambda.path), ")\n", sep=""))
	cat("\nLamba path : \n")
	print(x$lambda.path)
	
	cat("\nVariability :\n")
	print(signif(x$v, 4))
	
	cat("\nGraphs along the regularization path:\n")
	str(x$network)
}
#
plot.GMS <- function(x, fn="", th=1e-6, ...){
	cat("Plot optimal network\n")
	require(igraph)
	
	optNet <- x$network[[x$opt.index]]
	
	optNet[abs(optNet) > th] <- 1
	optNet[optNet!=1] <- 0
	diag(optNet) <- 0
	
	optG <- graph.adjacency(optNet,"undirected")
	#lout <- layout.fruchterman.reingold(optG)
	#lout <- layout.reingold.tilford(optG)
	
	if(is.null(colnames(optNet))){
		colnames(optNet) <- paste("V", 1:ncol(optNet), sep="")
	}

	V(optG)$label <- colnames(optNet)
	V(optG)$label.cex <- 35/nrow(optNet)
	#V(optG)$label.font <- 2
	V(optG)$label.color <- "#060606"
	#V(optG)$size <- 500/nrow(optNet)
	#V(optG)$label.cex = V(optG)$size*0.05
	
	V(optG)$size <- 900/nrow(optNet)
	V(optG)$label.cex = V(optG)$size*0.05
	V(optG)$label.font = 2
	
	V(optG)$frame.color <- NA
	V(optG)$shape <- "circle"
	V(optG)$color <- "#0099FF"
	
	E(optG)$width <- 50/nrow(optNet) * 2
	E(optG)$arrow.size <- 0
	E(optG)$curved <- 0.08
	E(optG)$color <- "#696A6A"
	
	if(fn != ""){
		pdf(fn, useDingbats = FALSE)
		#plot(optG, layout=layout.fruchterman.reingold(optG, niter=3000))
		plot(optG, layout=layout.kamada.kawai(optG, niter=1000))
		#plot(optG,layout=layout.drl)
		dev.off()
		cat(paste("Output file: ", fn, "\n",sep=""))
	}
	if(fn ==""){
		#plot(optG, layout=layout.fruchterman.reingold(optG, niter=3000))
		plot(optG, layout=layout.kamada.kawai(optG, niter=1000))
	}
}


# This will merge the net from each regression, for all lambdas
Merge.GraphXY <- function(network.X, network.Y, namesX, namesY, nlams, method="Both"){
	
	if (method == "Both"){
		results <- lapply(1:nlams, function(r){
						tmp1 <- network.X[,,r]
						rownames(tmp1) <- colnames(tmp1) <- c(namesX, namesY)
						
						tmp2 <- network.Y[,,r]
						rownames(tmp2) <- colnames(tmp2) <- c(namesY, namesX)
						
						tmp = tmp1
						tmp[namesX, colnames(tmp)] = tmp2[namesX, colnames(tmp)]
						return(tmp)})
	}
	
	if (method == "Right"){
		results <- lapply(1:nlams, function(r){
						tmp1 <- network.X[,,r]
						rownames(tmp1) <- colnames(tmp1) <- c(namesX, namesY)
						
						tmp2 <- network.Y[,,r]
						rownames(tmp2) <- colnames(tmp2) <- namesX
						
						tmp = tmp1
						tmp[namesX, namesX] = tmp2[namesX, namesX]
						return(tmp)})
	}
	
	if (method == "Left"){
		results <- lapply(1:nlams, function(r){
						tmp1 <- network.X[,,r]
						rownames(tmp1) <- colnames(tmp1) <- namesY
						
						tmp2 <- network.Y[,,r]
						rownames(tmp2) <- colnames(tmp2) <- c(namesY, namesX)
						
						tmp = tmp2
						tmp[namesY, namesY] = tmp1[namesY, namesY]
						return(tmp)})
	}
	
	return(results)
}

# Function
#	X is pxn
# 	Y is qxn
#	method 	- "Both": X <-> Y
#			- "Left": X <- Y
#			- "Right":  X -> Y
MGM.select <- function(X, Y, xlink="poisson", ylink="binomial", method="Both", N=100, beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL, parallel=T,nCpus=4, standardize=TRUE){
	 
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
#
#
# 	X is pxn = nvars x nobs 
#	Y default to NULL: if NULL, then do GM on X only. if not NULL, then Y is qxn
#	link is the method on Y
glmGeneric <- function(X, Y=NULL, link, lambda, parallel=F, nCpus = 4, standardize=TRUE){
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


# This function will create an empty "glmnet returned object" (*only certain elements*)
# Called in the case when glmnet returned an empty model and thrown an error
# Input: X is a n-obs by p-var matrix
glmnetEmpty <- function(X, lambda){
	fit = list()
	
	fit$a0 <- rep(0, length(lambda))
	fit$lambda <- lambda
	fit$df <- 0
	fit$dim <- c(ncol(X), length(lambda))
	
	fit$beta <- Matrix(0, ncol(X), length(lambda))
	rownames(fit$beta) <- colnames(X)
	colnames(fit$beta) <- paste("s", 0:(length(lambda)-1), sep="")
	
	return(fit)
}





