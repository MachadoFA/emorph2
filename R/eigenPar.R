#' A multivariate test of parallel evolution between multiple lineages
#'
#' Test for the presence of preferencial lines of evolution among multiple
#' independent lineages/replicates.
#'
#' @param X1 Matrix kxn for k number of traits and n lineages. A matrix of
#' multivariate averages of lineages before selection/
#' experimental manipulation OR a multivariate matrix of differences between
#' before and after selection.
#' @param X2 NULL or a kxn matrix of multivariate averages of lineages before
#' selection. If a value is provided, it has to be of the same dimentions as \code{X1}
#' @param null Character. Type of null hypothesis to be tested. If \code{null=
#' "random"}, tests if empirical vectors differs from totally random directions.
#' If null="drift", function tests if vectors differsn from what is expected
#' under drift. A \code{G} matrix has to be provided.
#' @param G Default to NULL. A kxk matrix of genetic covariances. Only used if
#' \code{null="drift"}.
#' @details The user can provide the function with either the lineages means before
#' and after selection/experimental manipulation (by using both \code{X1} and \code{X2}
#' respectively), or the vectors of differences between lineages.
#'
#' The function uses the within-lineage vectors of differences to calculates the
#'  matrix C of pairwise vector correlations
#'
#' \deqn{C=XX^t}
#'
#' with X being the matrix of within-lineage normalized vectors of divergence. The eigenvalues of C are a measure of how much divergence occured along the same directions.
#'
#' To evaluate if observed patterns cannot be explained by change alone,  the function employs a simulation procedure to generate a distribution of eigenvalues. The "random" eigenvalues are confronted agaist the observed ones to provide a p-value for each individual dimension. The function allows for the test of two null hypotheses. If \code{null="random"}, the null distribution is built by re-estimating C from completly random vectors, as suggested by De Lisle & Bolnick (2020). If \code{null="drift"}, vectors are simulated using the multivariate breeder's equation (Lande, 1979) as follows
#'
#' \deqn{\Delta x= G (t/Nef)}
#'
#' where \code{G} is the among traits genetic variance-covariance matrix, \code{t} is time and \code{Nef} is the effective population size. Because we are only interested in directions, both \code{t} and \code{Nef} can be ignored, but \code{G} has to be provided by the user.
#'
#' @return
#'  \describe{
#'      \item{eigenvalues}{ }
#'      \item{CI}{ }
#'      \item{A}{ }
#'      \item{vector correlation matrix}{ }
#'      }
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
