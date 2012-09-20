#### simulation code. 
#source('llgm.R')
source('PGMfunctions.R')

library(huge)

n = 100
p = 20
v = 0.1
u = 0.01
lambda = 0.3

gdata = huge.generator(n,d = p, graph="scale-free",v=v,u=u)

smatrix  = matrix(sample(c(1,-1), nrow(gdata$theta)*ncol(gdata$theta), replace =T), nrow = nrow(gdata$theta) )

simData = WPGMSim(n,p,R = 10, alpha = rep(0,p), Theta = .2*as.matrix(gdata$theta)*smatrix, maxit = 100 )


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


pdf('../Data/wpgm-hi.pdf')
myroc(as.matrix(gdata$theta),gpath[[1]],add=F)
myroc(as.matrix(gdata$theta),gpath[[2]],add=T,col="red")
myroc(as.matrix(gdata$theta),gpath[[3]]$icov,add=T,col="blue")
myroc(as.matrix(gdata$theta),gpath[[4]]$icov,add=T,col="green")
dev.off()

save.image(file="../Data/wpgm-hi.rdata")


