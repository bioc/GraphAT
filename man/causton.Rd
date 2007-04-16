\name{causton}
\alias{causton}
\docType{data}
\title{Yeast mRNA Expression Data}
\description{ This data set contains mRNA expression from a microarray experiment involving yeast grown under a variety of altered environments (e.g. acid, heat, sorbitol, etc.)  }
\usage{data(causton)}
\format{A matrix whose rows are the 6015 genes and whose columns are the 45 experimental conditions.}
\source{\url{http://web.wi.mit.edu/young/environment}}
\references{Causton, H. C. et al. 2001 ``Remodeling of Yeast Genome Expression in Response to Environmental Changes.'' Molecular Biology of the Cell \bold{12}, 323--337. 
}
\examples{
data(causton)

## Find the 3000 most variable genes, according to sd/mean:

varMeas<-function(vec) sd(vec)/mean(vec)
variability<-apply(causton,1,varMeas)

rks<-rank(variability)

causton3000<-causton[rks>length(rownames(causton))-3000,]
}
\keyword{datasets}
