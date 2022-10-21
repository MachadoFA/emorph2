#' Calculate interlandmark distances
#'
#' xxxxxxxx
#'
#' @param X Three dimentional array of form data.
#' @param dists Distances to be extracted
#' @param Names Landmark names
#' @details xxxxxxxxx
#' @return A dataframe of interlandmark distances per specimen.
#' @export
#' @references
#' @author Fabio Andrade Machado
#' @seealso \code{\link[ape]{pic}}
cdists<-function (X, dists, names=rownames(X)) {
	Dists<-matrix(NA,dim(X)[3],length(dists))

	dist.index<-t(matrix(unlist(strsplit(dists,"-")),2,length(dists)))

	l.index<- matrix(NA,length(names),2)
	l.names<-strsplit(names,"_")
	for (i in 1:length(names)) if (length(l.names[[i]])==1) l.index[i,1]<-l.names[[i]] 	else  l.index[i,]<-l.names[[i]]

	l.inv<-apply(dist.index,1,function(x) sum(l.index[,1]==x[1]|l.index[,1]==x[2]))

	for(i in 1:length(dists)){
		if(l.inv[i]==2){
			l1<-which(l.index[,1]==dist.index[i,1])
			l2<-which(l.index[,1]==dist.index[i,2])
			if(dist.index[i,1]==dist.index[i,2]){
			  Dists[,i]<- apply(X,3, function(x) sqrt(sum((x[l1[1],]-x[l2[2],])^2)))
			} else Dists[,i]<- apply(X,3, function(x) sqrt(sum((x[l1,]-x[l2,])^2)))
		}

		if(l.inv[i]==3){
			l.1<-which(l.index[,1]==dist.index[i,1])
			l.2<-which(l.index[,1]==dist.index[i,2])
			if(length(l.1)==2) {
				l2<-l.1
				l1<-l.2
			}
			if(length(l.1)==1) {
				l1<-l.1
				l2<-l.2
			}
			m1<-apply(X,3, function(x) sqrt(sum((x[l1,]-x[l2[1],])^2)))
			m2<-apply(X,3, function(x) sqrt(sum((x[l1,]-x[l2[2],])^2)))
			Dists[,i]<-apply(cbind(m1,m2),1,function(x) mean(x,na.rm=T))

		}
		if(l.inv[i]==4){
			l1<-which(l.index[,1]==dist.index[i,1])
			l2<-which(l.index[,1]==dist.index[i,2])

		m1<-apply(X,3, function(x) sqrt(sum((x[l1[l.index[l1,2]=="E"],]-x[l2[l.index[l2,2]=="E"],])^2)))
		m2<-apply(X,3, function(x) sqrt(sum((x[l1[l.index[l1,2]=="D"],]-x[l2[l.index[l2,2]=="D"],])^2)))
		Dists[,i]<-apply(cbind(m1,m2),1,function(x) mean(x,na.rm=T))
		}
	}
	return(Dists)
}
