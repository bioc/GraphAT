#function to make an adjacency matrix for a cluster graph
#nvec is a vector of the size of the clusters

makeClustM <- function(nvec){

	nClust<-length(nvec)

	c1mat<-1-diag(nvec[1])

	if(nClust>1){
	for (n in 2:nClust){
	
		c2mat<-1-diag(nvec[n])
	   	betweenmat<-matrix(0,nvec[n],sum(nvec[1:(n-1)]))
		c1mat<-cbind(rbind(c1mat,betweenmat),
				rbind(t(betweenmat),c2mat))
	}
	}
	return(c1mat)
}

