\name{GetDepthMat} 
\alias{GetDepthMat}


\title{Function to obtain the "depth" values for pairs of genes based
on each of the three GO networks (BP,CC,MF) }

\description{ This function takes as input a vector of YEAST probe
(gene) names as well as the specific ontology of interest, i.e. BP,CC,
or MF. The function returns a square matrix of size equal to the
number of genes in the input vector, where the i,jth element of the
output matrix corresponds to the maximum depth of all
annotations/terms shared in common between genes i and j with respect
to a specific GO graph. Depth of a specific annotation/term is defined
as the length of the shortest path between the term and the root node
of the ontology of interest.The diagonal elements of the output matrix
are set to 0. Note that if two genes share only the root node in
common, their depth is set to 1 (i.e.minimum of depth values for all
pairs of genes).}

\usage{GetDepthMat(genes.list, gocategory)} 

\arguments{
\item{genes.list}{Vector of YEAST gene (probe) names (labels).}
\item{gocategory}{Specific GO category of interest. Possible choices are: "BP","MF","CC"}
}
\value{A square matrix of size equal to the number of genes in input vector or gene (probe) names (labels).}
\examples{ 
genes.list <- c("YBR017C","YDL126C","YDL127W")
gocategory <- "CC"
depthmat <- GetDepthMat(genes.list,gocategory)
print(depthmat)
}
\author{Raji Balasubramanian \email{rbalasub@hsph.harvard.edu}}
\keyword{methods}
