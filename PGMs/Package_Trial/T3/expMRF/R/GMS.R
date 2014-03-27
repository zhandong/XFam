require(huge)
require(glmnet)
require(multicore)
require(igraph)
require(snowfall)

GMS <-
function(x, ...) UseMethod("GMS")

