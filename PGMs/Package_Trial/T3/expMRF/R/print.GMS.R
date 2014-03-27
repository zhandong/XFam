print.GMS <-
function(x, ...)
{
	cat("Call:\n")
	print(x$call)
	
	cat(paste("\nOptimal lambda index: ", x$opt.index, "\n", sep=""))
	
	#cat(paste("\nLamba path (length: ", length(x$lambda.path), ")\n", sep=""))
	cat("\nLamba path : \n")
	print(x$lambda.path)
	
	cat("\nVariability :\n")
	print(signif(x$v, 4))
	
	cat("\nGraphs along the regularization path:\n")
	str(x$network)
}
