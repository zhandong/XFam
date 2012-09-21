#### simulation code. 
#source('llgm.R')
source('PGMfunctions.R')

library(huge)

fooSim <-function(n = 200, p = 50, v =0.1, u = 0.01, lambda =0.3, theta = 0.1,fname="WPGM",graphType="scale-free")
{
	

# n = 200
# p = 50
# v = 0.1
# u = 0.01
# lambda = 0.3
# theta = 0.1

gdata = huge.generator(n,d = p, graph="scale-free",v=v,u=u)

smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace =T), nrow = nrow(gdata$theta) )

simData = WPGMSim(n,p,R = 10, alpha = rep(0,p), Theta = theta*as.matrix(gdata$theta)*smatrix, maxit = 100 )


######################
### WPGM
gpath =c()

# gpath[[1]] = WPGM.select(t(simData),R=max(simData),parallel=T,ncores = 20,nlams = 20,N=20)
# gpath[[2]] = LPGM.select(t(simData),nlams = 20, nCpus = 10, parallel=T)
# gpath[[3]]  = myglasso.select(t(simData),nlams=20)
# gpath[[4]]  = myglasso.select(t(huge.npn(simData)),nlams = 20)


gpath[[1]] = WPGM.network(t(simData),lmin = 0.01, R=max(simData),parallel=T,ncores = 20,nlams = 20)
gpath[[2]] = LPGM.network(t(simData),lmin=0.01,nlams = 20, nCpus = 10, parallel=T)
gpath[[3]]  = myglasso(t(simData),nlams=20)
gpath[[4]]  = myglasso(t(huge.npn(simData)),nlams = 20)


fname = paste(fname,"-n-",n,"-p-",p,"-theta-",theta,"-",graphType,sep="")

pdf(paste(fname,".pdf",sep=""))
myroc(as.matrix(gdata$theta),gpath[[1]],add=F)
myroc(as.matrix(gdata$theta),gpath[[2]],add=T,col="red")
myroc(as.matrix(gdata$theta),gpath[[3]]$icov,add=T,col="blue")
myroc(as.matrix(gdata$theta),gpath[[4]]$icov,add=T,col="green")
dev.off()

save.image(file=paste(fname,".rdata",sep=""))

}


for( theta in seq(0.1,0.2,by=0.02)){
	fooSim(n=200,p=50,theta = theta,graphType = "star");
}

