\name{clust2Mat}
\alias{clust2Mat}
\title{Function to compute adjacency matrix of cluster graph given a vector of cluster memberships}
\description{
  Given a list of cluster memberships, this function computes the adjacency matrix for the associated cluster graph.  The adjacency matrix is a matrix whose rows and columns are the nodes of the cluster graph and whose entries are 0's or 1's.  A 1 entry indicates that the corresponding nodes are connected, and a 0 indocates that they are not.  }
\usage{
clust2Mat(memb)
}
\arguments{
  \item{memb}{ A numeric vector, with each entry representing a node, the entry's value being the number of the cluster to which that node belongs.}} 
\details{
  Given a vector of cluster membership numbers, where the number of entries is the number of nodes n, the function computes an nxn ``adjacency matrix'' for the corresponding cluster graph.  The cluster graph is the graph in which two nodes are connected by an edge if and only if they are members of the same cluster.  The adjacency matrix for the graph has rows and columns representing the nodes, in the same order as the input vector.  The (i,j) entry is 1 if and only if node i and node j are in the same cluster.  Otherwise, the entry is 0.  By convention, diagonal entries are 0.  }
\value{An nxn adjacency matrix for the cluster graph, where n=length of cluster membership input vector memb.}
\author{Tom LaFramboise \email{tlaframb@hsph.harvard.edu}}
\seealso{\code{\link{makeClustM}}}
\examples{
memberships<-c(1,2,3,1,2,3,1,2,3,4)
clust2Mat(memberships)}
\keyword{data}
