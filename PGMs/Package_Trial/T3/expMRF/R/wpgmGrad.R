wpgmGrad <-
function(X,Y,R,beta)
{
	n = nrow(X); p = ncol(X);
	t2 = exp(matrix((0:R)%x%X%*%beta,n,R+1) - matrix(1,n,1)%*%t(matrix(log(factorial((0:R))),R+1,1)))
	denom = apply(t2,1,sum)
	t1 = array(t((0:R)%x%t(X)),c(n,p,R+1))
	t3 = array(matrix(1,p,1)%x%t2,c(n,p,R+1))
	num = apply(t1*t3,c(1,2),sum)
	return( -t(X)%*%Y + apply(num/(denom%*%matrix(1,1,p)),2,sum))
}
