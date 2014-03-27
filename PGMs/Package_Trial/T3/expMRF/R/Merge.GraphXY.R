Merge.GraphXY <-
function(network.X, network.Y, namesX, namesY, nlams, method="Both"){
	
	if (method == "Both"){
		results <- lapply(1:nlams, function(r){
						tmp1 <- network.X[,,r]
						rownames(tmp1) <- colnames(tmp1) <- c(namesX, namesY)
						
						tmp2 <- network.Y[,,r]
						rownames(tmp2) <- colnames(tmp2) <- c(namesY, namesX)
						
						tmp = tmp1
						tmp[namesX, colnames(tmp)] = tmp2[namesX, colnames(tmp)]
						return(tmp)})
	}
	
	if (method == "Right"){
		results <- lapply(1:nlams, function(r){
						tmp1 <- network.X[,,r]
						rownames(tmp1) <- colnames(tmp1) <- c(namesX, namesY)
						
						tmp2 <- network.Y[,,r]
						rownames(tmp2) <- colnames(tmp2) <- namesX
						
						tmp = tmp1
						tmp[namesX, namesX] = tmp2[namesX, namesX]
						return(tmp)})
	}
	
	if (method == "Left"){
		results <- lapply(1:nlams, function(r){
						tmp1 <- network.X[,,r]
						rownames(tmp1) <- colnames(tmp1) <- namesY
						
						tmp2 <- network.Y[,,r]
						rownames(tmp2) <- colnames(tmp2) <- c(namesY, namesX)
						
						tmp = tmp2
						tmp[namesY, namesY] = tmp1[namesY, namesY]
						return(tmp)})
	}
	
	return(results)
}
