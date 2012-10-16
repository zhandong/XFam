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


# Get some metrics on the network topology
GraphInfo <- function(tnet,  th=0.0001, type="adj"){
	if(type=="adj"){
			tnet[abs(tnet) > th] <- 1
			tnet[tnet!=1] <- 0
			diag(tnet) <- 0
			tnet <- graph.adjacency(tnet,"undirected")
	}
	
	subcomponent(tnet, V(tnet)[1], mode="all")

	# Get the density of the whole graph
	dens <- graph.density(tnet, loops=F)

	# Get the size (# of vertices) of the largest component
	comps <- sapply(1:vcount(tnet), function(i) length(subcomponent(tnet, V(tnet)[i], mode="all")))
	max.comp <- max(comps)
	max.comp.perc <- max.comp/vcount(tnet)

	# Get the density of the largest component
	max.comp.vs <- V(tnet)[comps==max.comp]
	max.comp.net <- induced.subgraph(tnet, max.comp.vs)
	max.comp.es <- ecount(max.comp.net)
	max.comp.dens <- graph.density(max.comp.net, loops=F)
	
	# Counts how many nodes are single
	sing.comp <- sum(comps==1)
	sing.comp.perc <- sing.comp/vcount(tnet)
	
	return(c(dens=dens, max.comp=max.comp, max.comp.perc=max.comp.perc, max.comp.es=max.comp.es, max.comp.dens=max.comp.dens, 
			sing.comp=sing.comp, sing.comp.perc=sing.comp.perc))
}












#
#
# The following part is generic - repeatable 
GeneralRun.glasso <- function(dat.exp.gd, fn.prefix){
	# Get optimal sparse net using glasso
	dat.gd.opt.glasso <- myglasso.select(dat.exp.gd, nlams=20)
	save(dat.gd.opt.glasso, file=paste(fn.prefix, ".saved", sep=""))

	dat.gd.opt.glasso.net <- dat.gd.opt.glasso$opt.icov
	dat.gd.opt.glasso.net <- as.matrix(dat.gd.opt.glasso.net)
	
	colnames(dat.gd.opt.glasso.net) <- rownames(dat.exp.gd)
	rownames(dat.gd.opt.glasso.net) <- rownames(dat.exp.gd)
	# save(dat.gd.opt.glasso.net, file=file=paste(fn.prefix, ".saved"))

	# Write the optimal net out
	GenGraphGml(dat.gd.opt.glasso.net,  fn=paste(fn.prefix, ".gml", sep=""), 0.00001, "adj")

	# Take a look at the network
	#dump <- apply(dat.gd.opt.glasso.net, 1, sum)
	#summary(dump)

	# Divide the optimal net to submodules, identified by colors
	dat.gd.opt.glasso.mod <- RandomWalk.Mod(dat.gd.opt.glasso.net, stp=10)
	GenGraphGml(dat.gd.opt.glasso.mod,  fn=paste(fn.prefix, ".mod.gml", sep=""), 0.00001, "igrpah")
}

#
#
# The following part is generic - repeatable 
GeneralRun <- function(dat.exp.gd, fn.prefix){
	# Get optimal sparse net using mb
	dat.gd.opt.mb <- mymb.select(dat.exp.gd, nlams=20)
	save(dat.gd.opt.mb, file=paste(fn.prefix, ".saved", sep=""))

	dat.gd.opt.mb.net <- dat.gd.opt.mb$refit
	colnames(dat.gd.opt.mb.net) <- rownames(dat.exp.gd)
	rownames(dat.gd.opt.mb.net) <- rownames(dat.exp.gd)
	# save(dat.gd.opt.mb.net, file=file=paste(fn.prefix, ".saved"))

	# Write the optimal net out
	GenGraphGml(dat.gd.opt.mb.net,  fn=paste(fn.prefix, ".gml", sep=""), 0.00001, "adj")

	# Take a look at the network
	#dump <- apply(dat.gd.opt.mb.net, 1, sum)
	#summary(dump)

	# Divide the optimal net to submodules, identified by colors
	dat.gd.opt.mb.mod <- RandomWalk.Mod(dat.gd.opt.mb.net, stp=10)
	GenGraphGml(dat.gd.opt.mb.mod,  fn=paste(fn.prefix, ".mod.gml", sep=""), 0.00001, "igrpah")
}

