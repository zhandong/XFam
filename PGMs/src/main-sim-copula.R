#### simulation code. 
source('llgm.R')
source('PGMfunctions.R')

library(huge)

n = 20
p = 10
v = 0.1
u = 0.01
lambda = 0.3

gdata = huge.generator(n,d = p, graph="scale-free",v=v,u=u)
simData = Copula.Norm.Pois(gdata$data,lambda = lambda)

######################
### WPGM
gpath =c()


gpath[[1]] = WPGM.select(t(simData),R=max(simData),parallel=T,ncores = 16,nlams = 20)
gpath[[2]] = LPGM.select(t(simData),nlams = 20, nCpus = 16, parallel=T)
gpath[[3]]  = myglasso.select(t(simData),nlams=20)
gpath[[4]]  = myglasso.select(t(huge.npn(simData)),nlams = 20)



pdf('copula-hi.pdf')
myroc(as.matrix(gdata$theta),gpath[[1]]$network,add=F)
myroc(as.matrix(gdata$theta),gpath[[2]]$network,add=F,col="red")
myroc(as.matrix(gdata$theta),gpath[[3]]$icov,add=T,col="blue")
myroc(as.matrix(gdata$theta),gpath[[4]]$icov,add=T,col="green")

dev.off()

save.image(file="copula-hi.rdata")


