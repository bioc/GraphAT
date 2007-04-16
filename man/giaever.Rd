\name{giaever}
\alias{giaever}
\docType{data}
\title{Yeast Gene-Knockout Fitness Data}
\description{ This data set contains fitness deficiency scores from gene knockout experiments involving yeast grown under a variety of altered environments (e.g. acid, heat, sorbitol, etc.)  }
\usage{data(giaever)}
\format{A matrix whose rows are the 5922 genes knocked out and whose columns are the 32 experimental conditions.}
\source{\url{http://gobi.lbl.gov/YeastFitnessData}}
\references{Giaever, G. et al. 2002 ``Functional profiling of the Saccharomyces cerevisiae genome.'' Nature \bold{418}, 387--391. 
}
\examples{
data(giaever)

## Find the 3000 most variable genes, according to sd/mean:

varMeas<-function(vec,na.rm=TRUE) 
     sd(vec,na.rm = na.rm)/mean(vec,na.rm = na.rm)

variability<-apply(giaever, 1, varMeas)

rks<-rank(variability)

giaever3000<-giaever[rks>length(rownames(giaever))-3000,]
}
\keyword{datasets}
