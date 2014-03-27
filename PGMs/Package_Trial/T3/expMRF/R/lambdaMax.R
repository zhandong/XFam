lambdaMax <-
function(X){
	tmp = t(X)%*%X
	return(max(	tmp[upper.tri(tmp)]))
}
