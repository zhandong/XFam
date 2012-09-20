#### load the miRNA 

source('funcs.R')
source('PGMfunctions.R')

load('../data/miRNAseq-breast.rdata')

#normalized by Bullard. 
ndata=t(normBullard(t(mydata)))

# filter out  the 25% quantile of the variance
index_var =apply(ndata,1,var)
index = index_var > quantile(index_var,0.50)

#normalize to poisson
Nmat = t(normG(t(ndata[index,]),quan=FALSE,nalphas=50))

#save.image('miRNAseq-data-brast-normalized.rdata')
#load("miRNAseq-data-brast-normalized.rdata")
Nmat = floor(Nmat)

#ghat = WPGM.network(Nmat,R=max(Nmat),nlams=5,F)

gpath = WPGM.select(Nmat)
save.image('WPGM-breast.rdata')

