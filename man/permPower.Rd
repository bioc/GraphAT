\name{permPower}
\alias{permPower}
\title{Function to compute estimated probability of detecting preferential
connection of intracluster nodes}
\description{
  This function simulates graphs from the alternative hypothesis
of preferential connection of intracluster nodes.  For each graph, it runs a node and edge permutation test.  The estimated ``power'' of each test is the proportion of graphs that the test rejects the null hypothesis of no preferential connection of intracluster edges.  }
\usage{
permPower(psi=1,clsizes, nedge, nhyper=100, nperms=1000)
}
\arguments{
  \item{psi}{The non-centrality parameter for the noncentral hypergeometric distribution used to simulate the graphs.}
\item{clsizes}{A vector of cluster sizes.}
\item{nedge}{The number of edges in each graph.}
\item{nhyper}{The number of noncentral hypergeometric graphs simulated to estimate "power".}
\item{nperms}{The number of permutations used for each run of the edge and node permutation tests.}
}
\details{
  The function first generates nhyper realizations of a noncentral hypergeometric(nedge,n,k,psi) random variable, where n is the number of node pairs and k is the number of intracluster node pairs.  For each realization x, a graph with n edges, x of which are intracluster, is generated.  The edge and node permutation tests (with nperms permutations each) are performed on each graph.  The estimated ``power'' of each test is the proportion of graphs for which the test rejects the null hypothesis of no preferential connection of intracluster nodes (at the 5\% level).  The 95\% confidence intervals for the power levels are also computed.  }
\value{A list with four components:
  \item{power.permedge}{Estimated ``power'' for edge permutation test.}
  \item{power.permnode}{Estimated ``power'' for node permutation test.}
\item{CI.permedge}{Vector giving 95\% confidence interval for edge permutation test power.}
\item{CI.permnode}{Vector giving 95\% confidence interval for node permutation test power.}
}
\author{Tom LaFramboise \email{tlaframb@hsph.harvard.edu}}
\seealso{\code{\link{permEdgesM2M}}, \code{\link{permNodesM2M}}, \code{\link{makeClustM}}}
\examples{
permPower(psi=5,clsizes=c(1,2,3,4),nedge=10,nhyper=100,nperms=100)
}
\keyword{htest}
