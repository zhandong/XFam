wpgmObj <-
function(X,Y,R,beta)
{
	n = nrow(X); p = ncol(X);
	t2 = exp(matrix((0:R)%x%X%*%beta,n,R+1) - matrix(1,n,1)%*%t(matrix(log(factorial((0:R))),R+1,1)))
	return( -t(Y)%*%X%*%beta + sum(log(apply(t2,1,sum))))
}
