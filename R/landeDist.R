#' Test for Lande's generalized genetic distance
#'
#' Test for the presence of preferencial lines of evolution among multiple
#' independent lineages/replicates.
#'
#' @param dz
#' @param G
#' @param Ne
#' @param CI
#' @param par
#' @param nsim
#' @param alpha
#' @param MonteCarlo
#' @param n
#' @param parallel
#' @details
#' @export
#' @references Lande, R. 1979. Quantitative genetic analysis of multivariate evolution, applied to brain: body size allometry. Evolution, 33(1), 402–416.
#' @author Fabio Andrade Machado
#' @importFrom mvtnorm rmvnorm
#'

landeDist.test <-
  function(dz, G, Ne, CI=0.95, par=FALSE, nsim=999, alpha=0.05, MonteCarlo=FALSE, n, parallel=FALSE){
  k<-ncol(dz)

  lD.observed<-landeDist(dz, G, Ne)

  if(par) {
    p <- pmin(2*(1-pchisq(lD.observed,df = k)),
              2*(pchisq(lD.observed,df = k)))
    thresh <- setNames(qchisq(c((1-CI)/2,(1+CI)/2),df = k),
                       c("2.5%","97.5%"))
  } else {
    `%dodo%` <- ifelse(parallel, `%dopar%`, `%do%`)
    lD.sim <- foreach(i=1:nsim,.combine = "c") %dodo% {
      if (MonteCarlo) {
        x<-rmvnorm(n,sigma=G)
        Gr<- var(x)
      } else {
        Gr <- G
      }

      z<-rmvnorm(1,sigma=G)
      landeDist(z,Gr,1)
    }

    p <- foreach(i=1:length(lD.observed), .combine=c) %do%
      min(mean(c(0.5,lD.sim<=lD.observed[i])),
          mean(c(0.5,lD.sim>=lD.observed[i])))*2
    thresh<-quantile(lD.sim, c((1-CI)/2,(1+CI)/2))

  }
  table<-data.frame(GGD=lD.observed, 'p-value'=p, 'consistent with'="drift")

  for(i in 1:nrow(table)){
    if(table[i,1]<thresh[1]) table[i,3]<-"stabilizing"
    if(table[i,1]>thresh[2]) table[i,3]<-"directional"
  }

  list(table=table, thresh=thresh)
}


#' Lande's generalized genetic distance
#'
#' Calculates Lande's generalized genetic distance
#'
#' @param dz
#' @param G
#' @param Ne
#' @details
#' @export
#' @references Lande, R. 1979. Quantitative genetic analysis of multivariate evolution, applied to brain: body size allometry. Evolution, 33(1), 402–416.
#' @author Fabio Andrade Machado

landeDist<-function(dz, G, Ne){
  mahalanobis(dz, center=FALSE, cov=G) * Ne
}

