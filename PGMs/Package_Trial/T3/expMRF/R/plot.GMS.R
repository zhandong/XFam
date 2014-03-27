plot.GMS <-
function(x, fn="", th=1e-6, ...){
	cat("Plot optimal network\n")
	require(igraph)
	
	optNet <- x$network[[x$opt.index]]
	
	optNet[abs(optNet) > th] <- 1
	optNet[optNet!=1] <- 0
	diag(optNet) <- 0
	
	optG <- graph.adjacency(optNet,"undirected")
	#lout <- layout.fruchterman.reingold(optG)
	#lout <- layout.reingold.tilford(optG)
	
	if(is.null(colnames(optNet))){
		colnames(optNet) <- paste("V", 1:ncol(optNet), sep="")
	}

	V(optG)$label <- colnames(optNet)
	V(optG)$label.cex <- 35/nrow(optNet)
	#V(optG)$label.font <- 2
	V(optG)$label.color <- "#060606"
	#V(optG)$size <- 500/nrow(optNet)
	#V(optG)$label.cex = V(optG)$size*0.05
	
	V(optG)$size <- 900/nrow(optNet)
	V(optG)$label.cex = V(optG)$size*0.05
	V(optG)$label.font = 2
	
	V(optG)$frame.color <- NA
	V(optG)$shape <- "circle"
	V(optG)$color <- "#0099FF"
	
	E(optG)$width <- 50/nrow(optNet) * 2
	E(optG)$arrow.size <- 0
	E(optG)$curved <- 0.08
	E(optG)$color <- "#696A6A"
	
	if(fn != ""){
		pdf(fn, useDingbats = FALSE)
		#plot(optG, layout=layout.fruchterman.reingold(optG, niter=3000))
		plot(optG, layout=layout.kamada.kawai(optG, niter=1000))
		#plot(optG,layout=layout.drl)
		dev.off()
		cat(paste("Output file: ", fn, "\n",sep=""))
	}
	if(fn ==""){
		#plot(optG, layout=layout.fruchterman.reingold(optG, niter=3000))
		plot(optG, layout=layout.kamada.kawai(optG, niter=1000))
	}
}
