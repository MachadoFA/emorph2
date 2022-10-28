#' @export
projectOnto<-function(A, target, tol=.Machine$double.eps){
  eig<-eigen(target)
  keep<-eig$values>tol
  if(isSymmetric(A)){
    X<-t(eig$vectors[,keep]) %*% A %*% eig$vectors[,keep]
  } else {
    X<-A %*% eig$vectors[,keep]
  }
  return(X)
}

#' @export
projectFrom<-function(A, target, tol=.Machine$double.eps){
  eig<-eigen(target)
  keep<-eig$values>tol
  if(isSymmetric(A)){
    X<-eig$vectors[,keep] %*% A %*% t(eig$vectors[,keep])
  } else {
    X<-A %*% t(eig$vectors[,keep])
  }
}
