# Simple create index function to use the ndIndex module
# 
# Author: tomf
###############################################################################
#' Create a new ndIndex
#' 
#' Create an index on a numeric matrix. Each row is a record, each column is a
#' dimension.
#' @param dataMatrix a numeric matrix object used to create the index.
#' @description Once created the findNeighbors property can be used to search
#'  this is a function that takes a point (with the same dimensionality as the
#'  index) and the number of points to return. 
#' @example 
#' #trivial example to show that distance and indexing works
#' index <- createIndex(matrix(c(0,0,3,4),ncol=2,byrow=T))
#' index$findNeighbors( c(0,3), 2, TRUE)
#' 
#' #now a million point index (takes about 2 seconds to build)
#' index <- createIndex(matrix(runif(2E6),ncol=2,byrow=T))
#' #find the 5 nearest points to (0.5, 0.5) (hopefully very, very fast!)
#' index$findNeighbors( c(0.5 , 0.5), 5, FALSE)
#' @export
createIndex<-function(dataMatrix){
	#make sure that the data is a matrix
	dataMatrix <- as.matrix(dataMatrix)
	
	#make sure that the data is numeric
	if(!is.numeric(dataMatrix)){
		stop("Indexes can only be created on numeric data.")
	}
	
	#create and return the index
	index <- new(ndIndex, dataMatrix)
	return(index)
}
