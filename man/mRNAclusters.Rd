\name{mRNAclusters}
\alias{mRNAclusters}
\docType{data}
\title{Yeast mRNA Expression Data Cluster Memberships}
\description{ This data set contains cluster membership for yeast genes clustered using mRNA expression from a microarray experiment in Causton, et al. Molecular Biology of the Cell (2001). The 3000 most variable genes were clustered using k-means with 30 clusters.}
\usage{data(mRNAclusters)}
\format{A data frame whose rows are the 3000 genes and whose two columns are gene name and cluster membership number.}
\source{\url{http://web.wi.mit.edu/young/environment}}
\references{Causton, H. C. et al. 2001 ``Remodeling of Yeast Genome Expression in Response to Environmental Changes.'' Molecular Biology of the Cell \bold{12}, 323--337. 
}
\examples{
data(mRNAclusters)

## Compute the adjacency matrix for the corresponding cluster graph:
mRNAMat<-clust2Mat(mRNAclusters[,2])
}
\keyword{datasets}
