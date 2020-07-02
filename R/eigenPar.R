#' A multivariate test of parallel evolution between multiple lineages
#'
#' Test for the presence of preferencial lines of evolution among multiple
#' independent lineages/replicates.
#'
#' @param X1 A matrix of multivariate averages of lineages before selection/
#' experimental manipulation OR a multivariate matrix of differences between
#' before and after selection.
#' @param X2 NULL or a matrix of multivariate averages of lineages before
#' selection. If a value is provided, it has to be of the same dimentions as X1
#' @param null Type of null hypothesis to be tested. If null="random", tests if
#' empirical vectors differs from totally random directions. If null="drift",
#' function tests if vectors differsn from what is expected under drift. A G
#' matrix has to be provided
#' @param G Default to NULL. A matrix of trait-by-trait genetic covariances.
#' Only used if null="drift".
#' @details
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom plyr aaply
#' @importFrom evolqg Normalize
#' @references De Lisle, S. P., & Bolnick, D. I. 2020. A multivariate view of parallel evolution. Evolution, 65, 1143–16. http://doi.org/10.1111/evo.14035
#' @author Fabio Andrade Machado
#' @seealso \code{\link[Morpho]{angleTest}}

eigenPar<-function(X1, X2=NULL, null=c("random", "drift"), sims=999, alpha=c(0.025,0.975), G=NULL, parallel=FALSE){
  if(!is.null(X2)){
    X1<-X1-X2
  }
  Xn<-apply(X1,1,Normalize)
  C <- t(Xn) %*% Xn
  eigenC<-eigen(C)

  null<-null[1]
  sim<- aaply(1:sims, 1, function(i){
    k <- dim(X1)
    if(null=="random")  M <- diag(k[2])
    if(null=="drift")   M <- G
    x <- rmvnorm(k[1], sigma = M)
    x <- apply(x,1,Normalize)
    c <- t(x)%*%x
    eigen(c)$values
  },.parallel = parallel)
  CI<-apply(sim,2,quantile,probs=alpha)

  A <- t(X1) %*% t(eigenC$vectors)
  colnames(A)<-paste("Axis",1:dim(A)[2])
  if(!is.null(rownames(X1))) rownames(A) <- rownames(X1)

  list('eigenvalues'=eigenC$values, CI=CI, A=A, 'vector correlation matrix'=C)
}
