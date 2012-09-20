# ------------------------------------------------------------------------------------------------------------------
#	Function definition
# 	Match gene expression and miR expression
GenemiR <- function(gtbl, mtbl){
	gpat <- data.frame(substr(colnames(gtbl)[-1], 1, 15))
	mpat <- data.frame(substr(colnames(mtbl)[-1], 1, 15))
	
	mexp <- data.frame(mpat, t(mtbl[,-1]))
	colnames(mexp) <- c("bcr", as.character(mtbl[,1]))
	
	gexp <- data.frame(gpat, t(gtbl[,-1]))
	colnames(gexp) <- c("bcr", as.character(gtbl[,1]))
	
	mgexp <- merge(mexp, gexp)
	return(mgexp)
}



# ------------------------------------------------------------------------------------------------------------------
# Load the mRNA and miRNA expression (Level 3) data
load("all.gene.saved")
load("all.mir.saved")

# ------------------------------------------------------------------------------------------------------------------
# Extract different groups of patients
#	tum.xx - patient samples with primary tumor
#	rec.xx - patient samples with recurrent tumor
#	nor.xx - normal samples

# mRNA expression
samplestypes <- sapply(1:length(names(all.gene)), function(i) unlist(strsplit(names(all.gene)[i], "-"))[4])
tum.all.gene <- cbind(all.gene[,1], all.gene[, grep("01", samplestypes)])
rec.all.gene <- cbind(all.gene[,1], all.gene[, grep("02", samplestypes)])
norm.all.gene <- cbind(all.gene[,1], all.gene[, grep("11", samplestypes)])
names(norm.all.gene)[1] <- names(all.gene)[1]

names(all.gene)[1] <- "gene.symbol"
names(tum.all.gene)[1] <- "gene.symbol"
names(rec.all.gene)[1] <- "gene.symbol"
names(norm.all.gene)[1] <- "gene.symbol"

tum.all.gene.bcr <- substr(colnames(tum.all.gene), 1, 12)

# miRNA expression
samplestypes <- sapply(1:length(names(all.mir)), function(i) unlist(strsplit(names(all.mir)[i], "-"))[4])
tum.all.mir <- cbind(all.mir[,1], all.mir[,grep("01", samplestypes)])
norm.all.mir <- cbind(all.mir[,1], all.mir[,grep("11", samplestypes)])
names(norm.all.mir)[1] <- names(all.mir)[1]

tum.all.mir.bcr <- substr(colnames(tum.all.mir), 1, 12)

# Get all gene names and miRNA symbols
genome <- as.character(tum.all.gene[,1])
MIR <- as.character(tum.all.mir[,1])


# Examples to combind gene and miRNA - this example done on patient samples with primary tumor
gene <- "LIN28"
gtbl <- tum.all.gene[tum.all.gene[,1] %in% gene,]
mtbl <- tum.all.mir

specgene.mir <- GenemiR(gtbl, mtbl)
head(specgene.mir[,1:5)
head(specgene.mir[,790:801])


