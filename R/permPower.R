
permPower<-function(psi=1,clsizes, nedge, nhyper=100, nperms=1000)
{

#function to  randomly select a subset of the edges of a graph
#takes an adjacency matrix called mat and number samp of edges to select
## Modified from Denise's permEdgesM2M

sampEdgesM2M<-function(mat, samp){

uniqE<-which(lower.tri(mat)==TRUE)
newE<-sample(uniqE,samp,replace=FALSE)
newmat<-matrix(0,dim(mat)[1],dim(mat)[2])
newmat[newE]<-1
newmat<-newmat+t(newmat)
return(newmat)
}


## First draw from noncentral hypergeometric distribution:
m1<-nedge
n1<-sum(choose(clsizes[clsizes>1],2))
n2<-choose(sum(clsizes),2)-n1
Xs<-rnoncenhypergeom(nhyper, n1, n2, m1, psi)

## Now make adjacency matrix for cluster graph and its complement
clustMat<-makeClustM(clsizes)
complMat<-(1-clustMat)-outer((1:sum(clsizes)),(1:sum(clsizes)),"==")

## Now a function to generate adjaceny matrix representing some X~noncenhyper
nCenMat<-function(x){
intraMat<-sampEdgesM2M(clustMat,x)
interMat<-sampEdgesM2M(complMat,m1-x)
return(intraMat+interMat)} 

### Now we define a function that returns a logical 2-vector indicating
### whether permedges and permnodes techniques, resp., reject null.
reject<-function(x){
randMat<-nCenMat(x)
cmpr.pedge<-function(num){
newM<-permEdgesM2M(randMat) ### Note that here we could use rhypergeom instead
return(x<=sum(newM*clustMat)/2)}
cmpr.pnode<-function(num){
newM<-permNodesM2M(randMat)
return(x<=sum(newM*clustMat)/2)}
return(c(mean(sapply(1:nperms,cmpr.pedge))<.05,
mean(sapply(1:nperms,cmpr.pnode))<.05))}

## The vector of whether the null was rejected. Top row for permuting edges
result<-sapply(Xs,reject)

## Return the approx. power and CI's
pws<-apply(result,1,mean)
return(list(power.permedge=pws[1],power.permnode=pws[2],
CI.permedge=binom.test(sum(result[1,]),dim(result)[2])$conf.int[1:2],
CI.permnode=binom.test(sum(result[2,]),dim(result)[2])$conf.int[1:2]))}

