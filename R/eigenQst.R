#' A multivariate Qst analysis
#'
#'
#' @param B A covariance matrix between population means.
#' @param G A genetic covariance matrix.
#'
#' @importFrom expm sqrtm
#' @author Fabio Andrade Machado
eigenQst<-function(B,G){
  mQst<-solve(sqrtm(B+2*G)) %*% B %*% solve(sqrtm(B+2*G))
  eQst<-eigen(mQst)
  globalQst<-mean(eQst$values)
  list('Global Qst'=globalQst,
       'Axes Qst'= eQst$values,
       'Axes of most selection'= eQst$vector)
}
