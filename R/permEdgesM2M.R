#function to permute the edges of a graph
#takes an adjacency matrix called mat

permEdgesM2M<-function(mat){

	numE<-sum(mat)/2
	uniqE<-which(lower.tri(mat)==TRUE)
	newE<-sample(uniqE,numE,replace=FALSE)
	newmat<-matrix(0,dim(mat)[1],dim(mat)[2])
	newmat[newE]<-1
	newmat<-newmat+t(newmat)
	return(newmat)
}
