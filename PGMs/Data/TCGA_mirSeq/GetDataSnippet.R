
setwd("/nfs/nfs_wooi/TCGA/OV/data/")
library(impute)
source("/nfs/nfs_wooi/TCGA/OV/src/preprocess-funcs.R")
source("/nfs/nfs_wooi/TCGA/OV/src/load.data.R")

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# ---
# miRNASeq:   tum.mirna.counts & tum.mirna.bcr
# miRNA array: tum.all.mir & tum.all.mir.bcr

dim(tum.mirna.counts)
dim(tum.all.mir)

# Identify bcr in both miRNA data
both.mirna.bcr <- intersect(tum.mirna.bcr, tum.all.mir.bcr)
length(both.mirna.bcr)

dat.gd <- tum.mirna.counts[,tum.mirna.bcr %in% both.mirna.bcr]

low.counts <- apply(dat.gd, 1, function(x) sum(x < 5)/ncol(dat.gd)) 
dat.gd.2 <- dat.gd[low.counts < 0.95, ]
dim(dat.gd.2)

ndata <- floor(t(normG(t(dat.gd.2),quan=TRUE,nalphas=20)))
filter2 <- apply(ndata, 1, var)
sum(filter2==0)
dat.gd.norm <- ndata[filter2 > 0, ]
dim(dat.gd.norm)


# Now, data (genes and samples) used to fit 
# original read counts: data.gd.2 
# read counts normalized to Poisson: dat.gd.norm (samed as that from table: /nfs/nfs_wooi/TCGA/OV/data/LPGM/G20.all/mirnaseq.gd.norm.G20.all.saved)





