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
