clust2Mat<-function(memb){
N<-length(memb)
return(as.numeric(outer(memb, memb, FUN="=="))-outer(1:N,1:N,"=="))}

