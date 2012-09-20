#### simulation code. 
source('llgm.R')
source('PGMfunctions.R')

n = 200
p = 50
alpha = rep(0,p)
mu = 0.5
sigma = sqrt(0.01)

#### randomly generate the graph adjacney and layout 
graphStr  = mp.generator(lambda=1, lambda.c=0.5, n, p, type="scale-free")
graphStr$layout  = myplot(graphStr$B)

Theta = matrix(rnorm(p*p,mu,sigma),nrow = p)
Theta = Theta * graphStr$B

############# WPGM sim

simData = WPGMSim(n,p,R=4,alpha,Theta,maxit = 100)
save.image("Scale-free-WPGM-hi-small.rdata")

load('Scale-free-WPGM-hi-small.rdata')
gpath = WPGM.select(t(simData),R=max(simData))


