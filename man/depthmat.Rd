\name{depthmat} 
\alias{depthmatBP}
\alias{depthmatCC}
\alias{depthmatMF}
\docType{data}

\title{Matrices of depth of association for pairs of YEAST genes with
respect to each of the BP,CC and MF ontologies of the GO database}

\description{ This matrix of depths is used to obtain the predictome
data in the paper. This is a symmetric matrix, where the i,j element
corresponds to is the maximum depth of all annotations shared by genes
i and j. Note that depth of a term in
a specific Gene Ontology (BP, CC, MF) is defined as the shortest path
between the term and the root node, where distance between nodes is
measured by the number of edges traversed. Row labels of the matrix
can be obtained by the row.names() function}

\usage{data(depthmatBP)} 

\format{Each of three matrices, namely depthmatBP.rda, depthmatCC.rda,
depthmatMF.rda is a symmetric matrix whose rows and columns correspond to
specific YEAST genes (see row labels using row.names()). The i,j entry
of each matrix refers to the maximum depth shared by genes i and j under each of the BP, CC and MF ontlogies respectively}

\source{\url{http://www.geneontology.org}}

\references{} 
\examples{ 
data(depthmatBP)
print(row.names(depthmatBP)[1:10])
}
\keyword{data}
