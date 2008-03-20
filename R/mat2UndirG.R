#V is list of node names
#mat is symmetric matrix indicating presence of edges

mat2UndirG<-function(V,mat){

	tmat<-which(mat==1,arr.ind=TRUE)
	numN<-length(V)
	numE<-dim(tmat)[1]

	rval <- vector("list", length = numN)
    	for (i in 1:numE) {
        	rval[[tmat[i, 1]]]$edges <- c(rval[[tmat[i, 1]]]$edges, 
            	tmat[i, 2])
        	ln <- length(rval[[tmat[i, 1]]]$edges)
        	rval[[tmat[i, 1]]]$weights <- c(rval[[tmat[i, 1]]]$weights, 
            	1)
        	names(rval[[tmat[i, 1]]]$weights)[ln] <- tmat[i, 2]
	}
    	names(rval) <- V
    	g1<-new("graphNEL", nodes = V, edgeL = rval)
	edgemode(g1)<-"undirected"
	return(g1)
}

