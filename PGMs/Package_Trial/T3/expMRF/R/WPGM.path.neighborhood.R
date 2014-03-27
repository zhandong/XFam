WPGM.path.neighborhood <-
function(X,Y,R,nlams,lmin=0.01,lambda=NULL,startb=0)
{
	n = nrow(X); p = ncol(X);
	
	if(is.null(lambda)){
		lmax = lambdaMax(t(X))
		lambda = exp(seq(log(lmax),log(lmin),l=nlams));	
	}
	
	if(nlams==1 & is.null(lambda)){
		lambda = lmax
	}
	thr = 1e-8; maxit = 1e6;
	Xt = cbind(t(t(rep(1,n))),X);
	if(sum(startb)==0){bhat = matrix(rnorm(p+1)/p,p+1,1)}else{bhat=startb}
	alphas = 0; Bmat = matrix(0,p,nlams);
	step = .1;
	for(i in 1:nlams){
		ind = 1; iter = 1;  
		while( thr<ind & iter<maxit){
			oldb = bhat; t = 1;
			grad = wpgmGrad(Xt,Y,R,oldb);
			oldobj = wpgmObj(Xt,Y,R,oldb);
			tmp = oldb - t*grad;
			bhat[1] = tmp[1]; 
			bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lambda[i]*t,max,0);
			newobj = wpgmObj(Xt,Y,R,bhat)
			
			while(newobj>9999999 | is.na(newobj) | is.na(newobj) ){
				t = t/p;
				tmp = oldb - t*grad;
				bhat[1] = tmp[1]; 
				bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lambda[i]*t,max,0);
				newobj = wpgmObj(Xt,Y,R,bhat)
			}
			
			while(newobj > oldobj - t(grad)%*%(oldb - bhat) + sum((oldb - bhat)^2)/(2*t)){
				t = t*step;
				tmp = oldb - t*grad;
				bhat[1] = tmp[1]; 
				bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lambda[i]*t,max,0);
				newobj = wpgmObj(Xt,Y,R,bhat)
			}

			iter = iter + 1;
			ind = sum((oldb - bhat)^2);
		}
		alphas[i] = bhat[1];
		Bmat[,i] = bhat[-1];
	}
	
	return(list(alpha=alphas,Bmat=Bmat,lambda=lambda))
}
