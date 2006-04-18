\name{getpvalue} 

\alias{getpvalue} 

\title{Function to obtain P values from the Edge permutation and Node
permutation tests respectively}

\description{The function takes as inputs two adjacency matrices. Let X denote the observed number of 
edges in common between the two adjacency matrices.  To
test the significance of the correlation between the two data sources,
the function performs N random edge permutations and random node
permutations respectively. For each permutation test, the function
outputs the proportion of N realizations that resulted in X edges or more at the intersection of the two datasources}

\usage{getpvalue(act.mat, nonact.mat, num.iterations = 1000)}

\arguments{ 

\item{act.mat}{Adjacency matrix corresponding to first data
source. That is, the i,j element of this matrix is 1 if data source
one specifies a functional link between genes i and j }

\item{nonact.mat}{Adjacency matrix corresponding to first data
source. That is, the i,j element of this matrix is 1 if data source
two specifies a functional link between genes i and j}

\item{num.iterations}{Number of realizations from random edge (node)
permutation to be obtained}
}
\details{We note that the first adjacency matrix, denoted act.mat is
the data source that is permutated with respect to edges or notes}

\value{A vector of length 2, where the first element is the P value
from Random Edge Permutation and the second element is the P value
from Random Node Permutation}

\author{Raji Balasubramanian \email{rbalasub@hsph.harvard.edu}}

\seealso{\code{\link{permEdgesM2M}}, \code{\link{permNodesM2M}}, \code{\link{makeClustM}}}
\examples{
act.mat <- matrix(0,3,3)
act.mat[2,1] <- 1 
act.mat[3,1] <- 1
nonact.mat  <- matrix(0,3,3)
nonact.mat[2,1] <- 1
nonact.mat[3,2] <- 1
p.val <- getpvalue(act.mat, nonact.mat, num.iterations = 100)
print(p.val)
}

\keyword{htest}
