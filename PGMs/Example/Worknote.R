setwd("C:/Users/yingwoow/Documents/GitHub/XFam/PGMs/Example")

# ---------------------------------
# Import necessary functions
# ---------------------------------
source("PGMfunctions.R")
source("preprocess-funcs.R")
source("postprocess-funcs.R")


# ---------------------------------
# Run on real data
# ---------------------------------

#####
# Load and preprocess data
load("mirSeqExp.saved")

# Filter out miRNA with no expression in samples
low.counts <- apply(dat.gd, 1, function(x) sum(x < 5)/ncol(dat.gd)) 
dat.gd.2 <- dat.gd[low.counts < 0.95, ]
dim(dat.gd.2)

# Normalize seq-data to Poisson distribution, with alpha = 20, 
# and remains miRNA that's variance != 0
ndata <- floor(t(normG(t(dat.gd.2),quan=TRUE,nalphas=20)))
filter2 <- apply(ndata, 1, var)
sum(filter2==0)
dat.gd.norm <- ndata[filter2 > 0, ]
dim(dat.gd.norm)


#####
# Start the analysis

# Set file variables to save output later
fn1 <- "LPGM/G20.all/mirnaseq.gd.norm.G20.all.saved"
fn2 <- "mirall.G20.all.path.saved"
fn3 <- "all.G20.all_path"
fn4 <- "all.G20.all.fitinfos.txt"

# Save normalized expression data
write.table(dat.gd.norm, file=fn1)

# Call LLGM model
# Short description on some paramters:
#		nlams 	: number of lambda for the regularization path, the larger the longer time to run
#		N 		: number of iterations for the stability selectio, the larger the longer time to run
#	Refer to function definition for other parameters
# mir.path.all = LPGM.select(dat.gd.norm, nlams=20, N=100, beta=0.05, nCpus=2, parallel=F)
mir.path.all = wooi.LPGM.select(dat.gd.norm, nlams=20, N=100, beta=0.05, nCpus=2, parallel=F)

# Save the results
save(mir.path.all, file=fn2)

# Generate grpahs and gather information
ind=20
its <- c(1:ind)
infos <- c()
for(i in its){
	temp.net <- as.matrix(mir.path.all$network[[i]])
	colnames(temp.net) <- rownames(dat.gd.norm)
	rownames(temp.net) <- rownames(dat.gd.norm)
	temp.info <- GraphInfo(temp.net)
	infos <- rbind(infos, temp.info)
	GenGraphGml(temp.net,  fn=paste(fn3, i, ".gml", sep=""), 0.00001, "adj")
}

select.infos <- data.frame(mir.path.all$v[its], mir.path.all$lambda.path[its])
rownames(select.infos) <- its
fit.infos <- data.frame(select.infos, infos)
write.table(fit.infos, file=fn4, sep="\t", quote=F)
