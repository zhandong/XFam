Copula.Norm.Pois <-
function(X,lambda)
{
	n = nrow(X); p = ncol(X);
	val = 0; dcuts = NULL; cnt = 0;
	while(val<max(.9999,1-2/(n*p))){
		val = ppois(cnt,lambda); cnt = cnt + 1;  dcuts = c(dcuts, val);
    }
	Y = matrix(0,n,p); oldval = min(X);
	for(i in 1:length(dcuts)) {
		val = quantile(X,dcuts[i]); Y[which(X<val & X>=oldval)] = i-1;
		oldval = val;
    }
	Y[X==max(X)] = max(Y) + 1;
  
	return(Y)
}
