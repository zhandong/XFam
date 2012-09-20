#################
# Poisson Graphical Models functions



##########
#PGM Gibbs Sampler - Besag's auto-model
#with negative relationships
#alpha a px1 vector
#Theta a pxp symmetric matrix (only off diags matter)
#note that only negative values in Theta are allowed
#n = sample size
#p = variable size
#maxit = iterations for Gibbs sampler

PGMSim = function(n,p,alpha,Theta,maxit)
{
  X = matrix(rpois(n*p,1),n,p);
  iter = 1;
  while(iter<maxit)
    {
      for(j in 1:p)
        {
          X[,j] = rpois(n,exp(alpha[j] + X[,-j]%*%Theta[-j,j]))
        }
      iter = iter + 1;
    }
  return(X)
}




##########
#Winsorized PGM Gibbs Sampler
#alpha a px1 vector
#Theta a pxp symmetric matrix (only off diags matter)
#n = sample size
#p = variable size
#R = truncating number
#maxit = iterations for Gibbs sampler

WPGMSim = function(n,p,R,alpha,Theta,maxit)
{
  X = matrix(rpois(n*p,1),n,p);
  iter = 1;
  while(iter<maxit)
    {
      for(j in 1:p)
        {
          num = exp( matrix(1,n,1)%*%t(alpha[j]*c(0:R)-log(factorial(c(0:R)))) + matrix(c(0:R)%x%X[,-j]%*%Theta[-j,j],n,R+1) );
          Pmat = num/matrix(apply(num,1,sum),n,R+1);
          X[,j] = apply(apply(Pmat,1,mymult)==1,2,which) - 1;
        }
      iter = iter + 1;
    }
  return(X)
}
  
#########
#internal
mymult = function(pvec){return(rmultinom(1,1,pvec))}



#################
#Copula transform from normal to Poisson
#takes
#X nxp matrix of Gaussians
#lambda = the Poisson mean for the transform
#spits out the copula transformed data matrix


Copula.Norm.Pois = function(X,lambda)
{
  n = nrow(X); p = ncol(X);
  val = 0; dcuts = NULL; cnt = 0;
  while(val<max(.9999,1-2/(n*p)))
    {
      val = ppois(cnt,lambda); cnt = cnt + 1;  dcuts = c(dcuts, val);
    }
  Y = matrix(0,n,p); oldval = min(X);
  for(i in 1:length(dcuts))
    {
      val = quantile(X,dcuts[i]); Y[which(X<val & X>=oldval)] = i-1;
      oldval = val;
    }
  Y[X==max(X)] = max(Y) + 1;
  return(Y)
}





###################
#function to compute paths for PGM neighborhood selection
#solves: min_{a,b} Y'*(X*b - a) + exp(a - X*b) + lam*|| b ||_{1}
#subject to b >= 0;
#note: we return negative b as this corresponds to LLGM 
#input:
#X  - nxp data matrix 
#Y - nx1 vector of responses
#nlams - number of lambdas for path (set nlams=1 to return form one
#value)
#lmax - maximum lambda value
#startb - defaul=0, otherwise a starting vector for [alpha beta']'
#output:
#alphas - 1 x nlams vector of intercepts
#Bmat - p x nlams sparse matrix of coefficients
#lams - the lambda values 


PGM.path.neighborhood = function(X,Y,nlams,lmax,startb=0)
{
  n = nrow(X); p = ncol(X);
  lams = exp(seq(log(lmax),log(.0001),l=nlams));
  if(nlams==1){lams = lmax};
  thr = 1e-6; maxit = 1e6;
  Xt = cbind(t(t(rep(1,n))),X);
  L = 1;
  alphas = 0; Bmat = matrix(0,p,nlams);
  if(sum(startb)==0){Bhat = matrix(runif(p+1),p+1,1)}else{Bhat=startb}
  for(i in 1:nlams){
    iter = 1; ind = 1;
    while(thr<ind & iter<maxit){
      oldb = Bhat;
      tmp = Bhat - (t(Xt)%*%Y - t(Xt)%*%exp(-Xt%*%Bhat))/L;
      Bhat = matrix(sapply(tmp - lams[i]/L,max,0),p+1,1); Bhat[1] = tmp[1];
      ind = sum((Bhat - oldb)^2)/sum((oldb^2));
      iter = iter + 1;
    }
    alphas[i] = Bhat[1];
    Bmat[,i] = -Bhat[2:(p+1),drop=FALSE];
  }
  return(list(alpha=alphas,Bmat=Bmat,lams=lams))
}




############
#WPGM neighborhood selection problem over a grid of lambdas
#X - nxp data
#Y - nx1 response
#nlams - number of regualrizaiton parameters
#lmax - max regularizaiton paramter
#startb - optional starting values
#Returns: alpha (intercept), beta (coefficient)

WPGM.path.neighborhood = function(X,Y,R,nlams,lmax,startb=0)
{
  n = nrow(X); p = ncol(X);
  lams = exp(seq(log(lmax),log(.0001),l=nlams));
  if(nlams==1){lams = lmax};
  thr = 1e-8; maxit = 1e6;
  Xt = cbind(t(t(rep(1,n))),X);
  if(sum(startb)==0){bhat = matrix(rnorm(p+1),p+1,1)}else{bhat=startb}
  alphas = 0; Bmat = matrix(0,p,nlams);
  step = .1;
  for(i in 1:nlams){
    ind = 1; iter = 1;  
    while( thr<ind & iter<maxit){
      oldb = bhat; t = 1;
      grad = wpgmGrad(Xt,Y,R,oldb);
      oldobj = wpgmObj(Xt,Y,R,oldb);
      tmp = oldb - t*grad;
      bhat[1] = tmp[1]; bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lams[i]*t,max,0);
      while(wpgmObj(Xt,Y,R,bhat) > oldobj - t(grad)%*%(oldb - bhat) + sum((oldb - bhat)^2)/(2*t))
        {
          t = t*step;
          tmp = oldb - t*grad;
          bhat[1] = tmp[1]; bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lams[i]*t,max,0);
        }
      iter = iter + 1;
      ind = sum((oldb - bhat)^2)/sum(oldb^2);
    }
    alphas[i] = bhat[1];
    Bmat[,i] = bhat[-1];
  }
  return(list(alpha=alphas,Bmat=Bmat,lams=lams))
}





############
#WPGM neighborhood selection problem
#X - nxp data
#Y - nx1 response
#lam - regualrizaiton parameter
#startb - optional starting values
#Returns: alpha (intercept), beta (coefficient)

WPGM.neighborhood = function(X,Y,R,lam,startb=0)
{
  n = nrow(X); p = ncol(X);
  thr = 1e-8; maxit = 1e6;
  Xt = cbind(t(t(rep(1,n))),X);
  if(sum(startb)==0){bhat = matrix(rnorm(p+1),p+1,1)}else{bhat=startb}
  step = .1; ind = 1; iter = 1;  
  while( thr<ind & iter<maxit){
    oldb = bhat; t = 1;
    grad = wpgmGrad(Xt,Y,R,oldb);
    oldobj = wpgmObj(Xt,Y,R,oldb);
    tmp = oldb - t*grad;
    bhat[1] = tmp[1]; bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lam*t,max,0);
    while(wpgmObj(Xt,Y,R,bhat) >  oldobj - t(grad)%*%(oldb - bhat) + sum((oldb - bhat)^2)/(2*t))
      {
        t = t*step;
        tmp = oldb - t*grad;
        bhat[1] = tmp[1]; bhat[-1] = sign(tmp[-1])*sapply(abs(tmp[-1]) - lam*t,max,0);
      }
    iter = iter + 1;
    ind = sum((oldb - bhat)^2)/sum(oldb^2);
  }
  return(list(alpha=bhat[1],beta=bhat[-1]))
}




#########
#gradient - internal
wpgmGrad = function(X,Y,R,beta)
{
  n = nrow(X); p = ncol(X);
  t2 = exp(matrix((0:R)%x%X%*%beta,n,R+1) - matrix(1,n,1)%*%t(matrix(log(factorial((0:R))),R+1,1)))
  denom = apply(t2,1,sum)
  t1 = array(t((0:R)%x%t(X)),c(n,p,R+1))
  t3 = array(matrix(1,p,1)%x%t2,c(n,p,R+1))
  num = apply(t1*t3,c(1,2),sum)
  return( -t(X)%*%Y + apply(num/(denom%*%matrix(1,1,p)),2,sum))
}


#############
#objective - internal
wpgmObj = function(X,Y,R,beta)
{
  n = nrow(X); p = ncol(X);
  t2 = exp(matrix((0:R)%x%X%*%beta,n,R+1) - matrix(1,n,1)%*%t(matrix(log(factorial((0:R))),R+1,1)))
  return( -t(Y)%*%X%*%beta + sum(log(apply(t2,1,sum))))
}


###################
#function to compute the maximum lambda
#X  - nxp data matrix 
lambdaMax <-function(X){
	
	tmp = t(X)%*%X
	return(max(	tmp[upper.tri(tmp)]))
}




################
#function to compute the poisson network
# X is pxn matrix
# nlams is the # of lambda
# parallel the algorithm if set T
WPGM.network = function(X,R,nlams,parallel=T){
	  lmax = lambdaMax(t(X))
	  lams = exp(seq(log(lmax),log(.0001),l=nlams));
	  ghat = c()
	  if(nlams>0){
	  		ghat = array(0,dim=c(nrow(X),nrow(X),length(lams)))
	   }
		wrapper <- function(i){
			
			print(i)
			
			fit = WPGM.path.neighborhood(t(X[-i,]),X[i,],R,nlams,lmax,0)
			fit$beta=as.matrix(fit$Bmat)
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
	ghat2=c()
	 if(parallel){

		library(multicore)
		ghat2=mclapply(1:nrow(X),wrapper)	
		}else{
		ghat2=lapply(1:nrow(X),wrapper)	
		}
		for(i in 1:nrow(X)){
			ghat[i,,]=ghat2[[i]]
		}
		
		return(ghat)
	
}	



################
#function to compute the poisson network
# X is pxn matrix
# nlams is the # of lambda
# parallel the algorithm if set T
PGM.network = function(X,nlams,parallel=T){
	  lmax = lambdaMax(t(X))
	  lams = exp(seq(log(lmax),log(.0001),l=nlams));
	  ghat = c()
	  if(nlams>0){
	  		ghat = array(0,dim=c(nrow(X),nrow(X),length(lams)))
	   }
		wrapper <- function(i){
			
			print(i)	
			fit = PGM.path.neighborhood(t(X[-i,]),X[i,],nlams,lmax,0)
			fit$beta=as.matrix(fit$Bmat)
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
	
	ghat2=c()
	if(parallel){

		library(multicore)
		ghat2=mclapply(1:nrow(X),wrapper)	
		}else{
		ghat2=lapply(1:nrow(X),wrapper)	
	}
	for(i in 1:nrow(X)){
			ghat[i,,]=ghat2[[i]]
	}
		
	return(ghat)
	
}	














