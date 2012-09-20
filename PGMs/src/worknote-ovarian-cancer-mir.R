##### generate the LPGM of ovarian miRNA network
source('funcs.R')
source('PGMfunctions.R')
source('llgm.R')

load('TCGA/ovarian.tumor.all.mir.rdata');

ovMir =tum.all.mir
rownames(ovMir)=ovMir[,1]
ovMir = ovMir[,-1]

#normalized by Bullard. 
ndata=t(normBullard(t(ovMir)))

# filter out  the 25% quantile of the variance
index_var =apply(ndata,1,var)
index = index_var > quantile(index_var,0.90)

#normalize to poisson
Nmat = t(normG(t(ndata[index,]),quan=FALSE,nalphas=50))

Nmat = floor(Nmat)

gpath = LPGM.select(Nmat,nlams=20,N=100,beta=0.05,nCpus=10)

optNet = gpath$network[[gpath$opt.index]]

rownames(optNet)=rownames(Nmat)
colnames(optNet)=rownames(Nmat)

library(igraph)

optNet[abs(optNet)>0.0001] = 1
optNet[optNet!=1]=0
diag(optNet)=0

write.graph(graph.adjacency(optNet,"undirected"),file="optNetOvarianMiRNA-90.gml",format = "gml")