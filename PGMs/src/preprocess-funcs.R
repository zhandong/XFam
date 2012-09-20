####################
#functions for miRNA preprocessing


normSeq = function(dat,ngams=5,nalphas=21,psi=5)
{
  gams = seq(.5,1,l=ngams)
  alphas = seq(.01,1, l=nalphas)
  GOFmat = matrix(0,length(gams),length(alphas))
  for(i in 1:length(gams)){
    for(j in 1:length(alphas)){
      qnum = quantile(rpois(1000,psi),gams[i])
      dis = apply(dat,1,quantile,gams[i])
      Nmatn = dat/(dis%o%rep(1,p)/qnum)
      betaj = apply(Nmatn,2,sum)
      dn = rep(1/n,n);
      GOFmat[i,j] = sum( (Nmatn^alphas[j] - dn%o%betaj )^2/(dn%o%betaj) )
    }
  }
  pmat = pchisq(GOFmat,(n-1)*(p-1),lower.tail=FALSE)
  indopt = which(pmat>.5)
  aopts = ceiling(indopt/ngams)
  alphaopt = alphas[max(aopts)]
  gopts = indopt%%ngams; gopts[gopts==0] = which.max(gams);
  gamopt = min(gams[gopts[aopts==max(aopts)]])
  
  qnum = quantile(rpois(1000,psi),gamopt)
  dis = apply(dat,1,quantile,gamopt)
  Nmatn = dat/(dis%o%rep(1,p)/qnum)
  x = Nmatn^alphaopt
  x = x/(apply(x,1,mean)%o%rep(1,p)/psi)    
  return(x)
}


normBullard = function(dat,psi=5)
{
  qnum = quantile(rpois(1000,psi),.75)
  dis = apply(dat,1,quantile,.75)
  p=ncol(dat)
  x = dat/(dis%o%rep(1,p)/qnum)
  return(x)
}

normCom = function(dat,nalphas=21,psi=5)
{
  qnum = quantile(rpois(1000,psi),.75)
  dis = apply(dat,1,quantile,.75)
  p=ncol(dat) #### Please check! 
  Nmat = dat/(dis%o%rep(1,p)/qnum)
  alphas = seq(.01,1, l=nalphas)
  betaj = apply(Nmat,2,sum)
  dis = apply(Nmat,1,sum)/sum(betaj)
  GOFmat = NULL
  for(j in 1:length(alphas)){
    GOFmat[j] = sum( (Nmat^alphas[j] - dis%o%betaj )^2/(dis%o%betaj) )    
  }
  alphaopt = alphas[which.min(GOFmat)]
  x = Nmat^alphaopt
  return(x)
}



normLi = function(dat,nalphas=21,psi=10)
{
  alphas = seq(.01,1, l=nalphas)
  betaj = apply(dat,2,sum)
  dis = apply(dat,1,sum)/sum(betaj)
  GOFmat = NULL
  for(j in 1:length(alphas)){
    GOFmat[j] = sum( (dat^alphas[j] - dis%o%betaj )^2/(dis%o%betaj) )    
  }
  alphaopt = alphas[which.min(GOFmat)]
  Nmat = dat^alphaopt
  x = Nmat/(apply(Nmat,1,mean)%o%rep(1,p)/psi)    
  return(x)
}



normG = function(dat,quan=TRUE,nalphas=50)
{  
  p = ncol(dat); n = nrow(dat)
  if(quan){
    dis = apply(dat,1,quantile,.75)
    qnum = mean(dis)
    dat = dat/(dis%o%rep(1,p)/qnum)    
  }
  alphas = seq(.01,1, l=nalphas)
  GOFmat = NULL
  for(j in 1:length(alphas)){
    GOFmat[j] = sum(abs(apply(dat^alphas[j],1,mean) - apply(dat^alphas[j],1,var))/apply(dat^alphas[j],1,mean))
  }
  alphaopt = alphas[which.min(GOFmat)]
  x = dat^alphaopt
  return(x)
}
