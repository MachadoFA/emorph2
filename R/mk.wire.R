#' @export
mk.wire<-function (X, dists, names=rownames(X)) {
  dist.index<-t(matrix(unlist(strsplit(dists,"-")),2,length(dists)))

  l.index<- matrix(NA,length(names),2)
  l.names<-strsplit(names,"_")
  for (i in 1:length(names)) if (length(l.names[[i]])==1) l.index[i,1]<-l.names[[i]] 	else  l.index[i,]<-l.names[[i]]

  l.inv<-apply(dist.index,1,function(x) sum(l.index[,1]==x[1]|l.index[,1]==x[2]))
  wire<-NULL
  # dinms<-NULL
  for(i in 1:length(dists)){
    if(l.inv[i]==2){
      l1<-which(l.index[,1]==dist.index[i,1])
      l2<-which(l.index[,1]==dist.index[i,2])
      wire<-rbind(wire,c(l1,l2))
      #dinms<-c(dinms,paste(l.index[l1,1],l.index[l2,1],sep="-"))
      }

    if(l.inv[i]==3){
      l.1<-which(l.index[,1]==dist.index[i,1])
      l.2<-which(l.index[,1]==dist.index[i,2])
      if(length(l.1)==2) {
        l1<-l.1
        l2<-c(l.2,l.2)
        #dinms<-c(dinms,paste(paste(l.index[l.1,1],l.index[l.1,2], sep="_"), l.index[l.2,],sep="-"))
        }
      if(length(l.1)==1) {
        l1<-c(l.1,l.1)
        l2<-l.2
        # dinms<-c(dinms,paste(l.index[l.1,],paste(l.index[l.2,1],l.index[l.2,2], sep="_"),sep="-"))
        }
      wire<-rbind(wire,cbind(l1,l2))
      }

    if(l.inv[i]==4){
      l1<-which(l.index[,1]==dist.index[i,1])
      l2<-which(l.index[,1]==dist.index[i,2])
      # dinms<-c(dinms,paste(paste(l.index[l.1,1],l.index[l.1,2], sep="_"),
      #                      paste(l.index[l.2,1],l.index[l.2,2], sep="_"),sep="-"))
      wire<-rbind(wire,cbind(l1,l2))
      }
  }
  # dinms<-sub("_NA","",dinms)
  # rownames(wire)<-dinms
  return(wire)
}
