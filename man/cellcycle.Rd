\name{cellcycle}
\alias{ccCM}
\alias{nNamescc}
\docType{data}
\title{Cell-Cycle Cluster Matrix}
\description{
An adjacency matrix in which 
}
\usage{data(ccCM)}
\format{
\code{ccCM} is a symmetric matrix with 2885 columns and 2885 rows.

\code{nNamescc} is a vector of 2885 gene names.
}
\details{
Cho, et al. discuss the k means clustering of 2885 Saccharomyces genes into 30 clusters with measurements taken over two synchronized cell cycles.  \code{nNamescc} is a vector of the 2885 gene names.  \code{ccCM} is an adjacency matrix in which a "1" in the ith row and jth column indicates that gene i and gene j belong to the same cluster.  All other entries are 0.  These data are integrated with phenotypic data and GO data in Balasubramanian, et al (2004).  

}
\source{
Balasubramanian R, LaFramboise T, Scholtens D, Gentleman R. (2004) A graph theoretic approach to integromics - integrating disparate sources of functional genomics data
}
\references{
Cho, et al. (1998) A genome-wide transcriptional analysis of the mitotic cell cycle.  Molecular Cell, 2, 65-73.

Tavazoie, et al. (1999) Systematic determination of genetic network architecture. Nature Genetics, 22, 281-285.
}
\examples{
data(ccCM)

}
\keyword{datasets}
