#' Lande's generalized genetic distance
#'
#' Test for the presence of preferencial lines of evolution among multiple
#' independent lineages/replicates.
#'
#' @param dz
#' @param G
#' @param Ne
#' @details
#' @export
#' @references Lande, R. 1979. Quantitative genetic analysis of multivariate evolution, applied to brain: body size allometry. Evolution, 33(1), 402–416.
#' @author Fabio Andrade Machado
landeDist<-function(dz, G, Ne){
  diag(dz %*% solve(M) %*% t(dz))*Ne
}

landeDist.test<-function(dz, G, Ne, par=FALSE, nsim=999, MonteCarlo=FALSE,parallel=F){
  k<-nrow(dz)

  lD<-landeDist(dz, G, Ne)
  if(par) {
    p<-pchisq(lD,df = k, lower.tail=TRUE)+pchisq(lD,df = k, lower.tail=TRUE)/2
  } else {
    for(i in 1:nsim){
      rmvnorm(dim(dat)[1],sigma=G)
      original<-landeDist(change,G,1)
      x<-mvtnorm::rmvnorm(dim(dat)[1],sigma=G)
      Gr<- var(x)
      resampled<-landeDist(change,Gr,1)
      Ge<-ExtendMatrix(Gr,ret.dim = 10)$ExtMat
      extended<-landeDist(change,Ge,1)
      data.frame(original, resampled, extended)
    }
  }


}
