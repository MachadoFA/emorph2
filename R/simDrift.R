#' Simulation of drift in a metapopulation..
#'
#' Simulates evolution by genetic drift on a metapopulation with varying levels
#' of geneflow.
#
#' @param G matrix kxk for k number of traits.
#' @param theta matrix of among-population coancestry coefficients.
#'
#' @return
#' A nxkxnsims cubic array of simulated population averages.
#'
#' @importFrom mvtnorm rmvnorm
#' @export

simDrift <-function(G, theta, nsims=1, mu=matrix(0,nrow(theta),nrow(G))){
  if(nrow(G)!=ncol(G)) stop("G is not square")
  if(!isSymmetric(G)) stop("G is not symmetric")
  if(nrow(theta)!=ncol(theta)) stop("theta is not square")
  if(!isSymmetric(theta)) stop("theta is not symmetric")

  k<-nrow(G)
  n<-nrow(theta)

  Sigma <- 2*G %x% theta

  x<-rmvnorm(nsims, sigma = Sigma)
  array(x, c(n,k,nsims))
}


