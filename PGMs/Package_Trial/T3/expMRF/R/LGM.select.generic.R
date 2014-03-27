LGM.select.generic <-
function(X, method="LPGM", link="poisson", N=100,beta=0.05, lmin = 0.01, nlams=20, lambda.path=NULL ,parallel=TRUE,nCpus=4){
	require('huge')
	require('glmnet')
	
	if(is.null(lambda.path) ){
		#lmax = myglmnet.max(X)
		lmax = myglmnet.max(X, link=link)
 		lambda.path = exp(seq(log(lmax),log(lmin),l=nlams));
	}
	
	if(parallel == T){
		b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
		ghat=list()
		ghat.path=list()
		ghat.path$path=vector("list",length(lambda.path))
		v=c()

		for(i in 1:N){
			cat(paste(method, ": Conducting sampling ... in progress: ", floor(100*(i/N)), "%", collapse=""),"\r")
			flush.console()
			
			glmpois.good <- 1
			
			while(glmpois.good){
				# Make sure sample with no gene with all zero values
				good <- 1
				while(good){
					index = sample(1:ncol(X),b,replace=F)
					#-- if(sum(apply(X[,index], 1, sum)==0)==0){
					if(sum(apply(X[,index], 1, function(x) length(unique(x))==1))==0){
						good <- 0
					}
				}
				
				tryCatch(
						{
							#ghat.path$raw= glmpois(X[,index],lambda=lambda.path,parallel=T,nCpus=nCpus)
							ghat.path$raw= glmGeneric(X[,index], NULL, link=link, lambda=lambda.path,parallel=T,nCpus=nCpus)	
							glmpois.good <- 0
						},
						error = function(e) {
							cat("glmnet returns empty model. Try again.")
						}
				)
			}
			
			for(j in 1:length(lambda.path)){
				tmp=ghat.path$raw[,,j]
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
		
		v=cummax(v)
		ghat$v=v
		ghat$lambda.path = lambda.path
		ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]	
			
		#ghat$network = glmpois(X,lambda=lambda.path,parallel=T,nCpus=nCpus)
		ghat$network = glmGeneric(X, NULL, link=link, lambda=lambda.path,parallel=T,nCpus=nCpus)
		ghat$network =lapply(1:nlams,function(r){return(ghat$network[,,r])})
		ghat$opt.index = which(v==max(v[v<beta]))
		ghat$call <- match.call()
		cat(paste("\n", method, " Completed.", "\n", sep=""))
		class(ghat) <- "GMS"
		return(ghat)
	}
	
	
	if(parallel == F){
		b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
		ghat=list()
		v=c()
		
		for( j in 1:length(lambda.path)){
			#cat ("j=", j, " \t")
			cat(paste(method, ": Conducting sampling ... in progress: ", floor(100*(j/length(lambda.path))), "%", collapse=""),"\r")
			flush.console()
			D=matrix(0,nrow=nrow(X),ncol=nrow(X))
			
			for(i in 1:N){
				#cat("\n i=", i, "\t")
				glmpois.good <- 1
				
				while(glmpois.good){
					# Make sure sample with no gene with all zero values
					good <- 1
					while(good){
						index = sample(1:ncol(X),b,replace=F)
						if(sum(apply(X[,index], 1, function(x) length(unique(x))==1))==0){
							good <- 0
						}
					}
					
					tryCatch(
						{
							#tmp=glmpois(X[,index],lambda=lambda.path[j],parallel=F)
							tmp=glmGeneric(X[,index],NULL, link=link, lambda=lambda.path[j],parallel=F)
							glmpois.good <- 0
						},
						error = function(e) {
							cat("glmnet returns empty model. Try again.\n")
						}
					)
				} 
				
				tmp[abs(tmp)<1e-06]=0
				tmp[abs(tmp)>1e-06]=1
				D=D+tmp
			}
			
			D=D/N
			D=2*D*(1-D)
			v=c(v,mean(D[upper.tri(D)]))			
		}
		
		v=cummax(v)
		ghat$v=v
		ghat$lambda.path = lambda.path
		ghat$opt.lambda = lambda.path[which(v==max(v[v<beta]))]	
		#ghat$network = glmpois(X,lambda=lambda.path,parallel=parallel,nCpus=nCpus)
		ghat$network = glmGeneric(X,NULL, link=link, lambda=lambda.path,parallel=parallel,nCpus=nCpus)
		ghat$network =lapply(1:nlams,function(r){return(ghat$network[,,r])})
		ghat$opt.index = which(v==max(v[v<beta]))
		ghat$call <- match.call()
		cat(paste("\n", method, " Completed.", "\n", sep=""))
		
		class(ghat) <- "GMS"
		return(ghat)
	}
		
}
