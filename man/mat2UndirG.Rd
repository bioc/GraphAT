\name{mat2UndirG}
\alias{mat2UndirG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Change graph representations}
\description{
A function to turn an adjacency matrix for a graph into a graphNEL object.
}
\usage{
mat2UndirG(V, mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{A vector of node names}
  \item{mat}{A square symmetric matrix indicating the presence of edges}
}
\details{
\code{mat} is a square matrix with rows and columns corresponding to nodes in the graph.  Entries of "0" indicate the lack of an edge.  Since this is making an undirected graph, \code{mat} must be symmetric.

}
\value{
A graphNEL object.
}
\references{Balasubramanian, et al. (2004) A graph theoretic approach to integromics - integrating disparate sources of functional genomics data.}
\author{Denise Scholtens}

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{
library(graph)
a <- matrix(c(0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0),ncol=4)
ag <- mat2UndirG(V=letters[1:4],mat=a)
}
\keyword{graphs}

