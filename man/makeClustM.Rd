\name{makeClustM}
\alias{makeClustM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Make an adjacency matrix for a cluster graph}
\description{
This function takes a vector of cluster sizes and returns an adjacency matrix for a graph in which edges connect nodes if they are members of the same cluster.
}
\usage{
makeClustM(nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nvec}{A vector of cluster sizes}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
A square adjacency matrix with the number of rows and columns equal to the sum of nvec.  An entry of "1" in the ith row and jth column indicates that node i and node j are members of the same cluster.  All other entries are "0".
}
\references{Balasubramanian, et al. (2004) A graph theoretic approach to integromics - integrating disparate sources of functional genomics data.}
\author{Denise Scholtens }
%\note{ ~~further notes~~ }

 %~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{clust2Mat}}}
\examples{
a <- makeClustM(c(2,3,4))
}
\keyword{graphs}

