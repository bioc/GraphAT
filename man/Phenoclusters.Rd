\name{Phenoclusters}
\alias{Phenoclusters}
\docType{data}
\title{Yeast Gene-Knockout Fitness Data Cluster Memberships}
\description{ This data set contains cluster memberships for yeast genes clustered using fitness deficiency scores from gene knockout experiments from Giaever et al. Nature (2002). The 3000 most variable genes were clustered using k-means with 30 clusters}
\usage{data(Phenoclusters)}
\format{A matrix whose rows are the 3000 genes and whose two columns are gene name and cluster membership number.}
\source{\url{http://gobi.lbl.gov/YeastFitnessData}}
\references{Giaever, G. et al. 2002 ``Functional profiling of the Saccharomyces cerevisiae genome.'' Nature \bold{418}, 387--391. 
}
\examples{
data(Phenoclusters)

## Compute the adjacency matrix for the corresponding cluster graph:
phenoMat<-clust2Mat(Phenoclusters[,2])

}
\keyword{datasets}
