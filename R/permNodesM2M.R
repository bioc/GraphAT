#permutes rows and columns in matrix
#equivalent to node permutation in a graph

permNodesM2M<-function(mat){
	
	ord<-sample(dim(mat)[1],dim(mat)[1],replace=FALSE)
	newmat<-mat[ord,ord]
	return(newmat)
}

