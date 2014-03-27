WPGM.neighborhood <-
function(X,Y,R,lam,startb=0)
{
	n = nrow(X); p = ncol(X);
	thr = 1e-8; maxit = 1e6;
	Xt = cbind(t(t(rep(1,n))),X);
	if(sum(startb)==0){
		bhat = matrix(rnorm(p+1)*.01,p+1,1)
	}
	else{
		bhat=startb
	}
	
	step = .1; ind = 1; iter = 1;  
	while( thr<ind & iter<maxit){
		oldb = bhat; t = 1;
		grad = wpgmGrad(Xt,Y,R,oldb);
		oldobj = wpgmObj(Xt,Y,R,oldb);
		tmp = oldb - t*grad;
		bhat[1] = tmp[1]; bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lam*t,max,0);
		while(wpgmObj(Xt,Y,R,bhat) >  oldobj - t(grad)%*%(oldb - bhat) + sum((oldb - bhat)^2)/(2*t)){
			t = t*step;
			tmp = oldb - t*grad;
			bhat[1] = tmp[1]; 
			bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lam*t,max,0);
		}
		iter = iter + 1;
		ind = sum((oldb - bhat)^2)/sum(oldb^2);
	}
	return(list(alpha=bhat[1],beta=bhat[-1]))
}
