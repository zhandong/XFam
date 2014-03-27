WPGM.select <-
function(X,R=max(X),N=100,beta=0.05,lmin=0.0001, nlams=20, lambda.path=NULL, parallel=FALSE, ncores = 4){
	require('huge')
	require('glmnet')
	
	if(is.null(lambda.path) ){
		lmax = lambdaMax(t(X))
 		lambda.path = exp(seq(log(lmax),log(lmin),l=nlams));
	}
	b = min(c(10*sqrt(ncol(X)), 0.8*ncol(X))) 
	ghat=list()
	ghat.path=list()
	ghat.path$path=vector("list",length(lambda.path))
	v=c()
	
	for(i in 1:N){
		cat(paste("WPGM: Conducting sampling ... in progress: ", floor(100*(i/N)), "%", collapse=""),"\r")
		flush.console()
		index = sample(1:ncol(X),b,replace=F)
		#tmp=glmpois(X[,index],lambda.path[j],parallel=parallel,warmStart=warmStart,nCpus=nCpus)
		#ghat.path$raw = WPGM.network(X[,index],R,nlams=length(lambda.path),parallel=parallel ,ncores = ncores)
		ghat.path$raw = WPGM.network(X[,index],R,nlams=length(lambda.path),lambda=lambda.path,parallel=parallel ,ncores = ncores)
				
		for(j in 1:length(lambda.path)){
			tmp=ghat.path$raw[[j]]
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
	ghat$network= WPGM.network(X,R,nlams=length(lambda.path),lambda=lambda.path,parallel=T)
	ghat$opt.index = which(v==max(v[v<beta]))

	cat("\nWPGM Completed. \n")

	return(ghat)
}
