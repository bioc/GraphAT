\name{perms}
\alias{permEdgesM2M}
\alias{permNodesM2M}
\title{Randomly permute edges or node labels in a graph}
\description{
Given an adjacency matrix for a graph, \code{permEdgesM2M} will return an adjacency matrix after an Erdos-Renyi random permutation of the edges in the graph.  \code{perNodesM2M} will return an adjacency matrix for a graph with identical structure, but with the node labels permuted.
}
\usage{
permEdgesM2M(mat)
permNodesM2M(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{A square adjacency matrix for a graph.}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
A square adjacency matrix for the new graph, subject to a random permutation of the edges or nodes.
}
\references{Balasubramanian, et al. (2004) A graph theoretic approach to integromics - integrating disparate sources of functional genomics data.}
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{permPower}}}
\examples{
g <- matrix(c(0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0),nrow=4)

g1 <- permEdgesM2M(g)
g2 <- permNodesM2M(g)

}
\keyword{graphs}

