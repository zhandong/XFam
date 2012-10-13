library(RColorBrewer)

####################
#functions for generating outputs on the optimal networks generated

# Generate gml file of the network 
# optNet is pxp adjacency matrix by default
# if "type" is specified as "igraph", then gml file is generated directly
# if "type" is specified otherwise, the adjacency matrix will obtain, then generate a igraph
GenGraphGml <- function(optNet, fn, th=0.0001, type="adj"){
	if(type!="igrpah"){
		if(type=="adj"){
			optNet[abs(optNet) > th] <- 1
			optNet[optNet!=1] <- 0
			diag(optNet) <- 0
			optNet <- graph.adjacency(optNet,"undirected")
		}
	}
	write.graph(optNet, file=fn , format="gml")
}

# Use random walk to modularize a networks into 
# input optNet is pxp adjacency matrix
RandomWalk.Mod <- function(optNet, stp=10){
	isol <- apply(optNet, 1, sum)
	optNet <- optNet[isol>0, isol>0]

	opt.g <- graph.adjacency(optNet, "undirected", weighted=T)
	rand.mod.net <- walktrap.community(opt.g,steps=stp)

	opt.g <- set.vertex.attribute(opt.g, "membership", value=rand.mod.net$membership)
	mycol <- rainbow(length(unique(rand.mod.net$membership)))
	mycol <- apply(col2rgb(mycol),2, function(x) paste(x, collapse=","))
	cols <- mycol[rand.mod.net$membership]
	opt.g <- set.vertex.attribute(opt.g, "color", value=cols)
	return(opt.g)
}
