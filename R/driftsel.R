#' Multivariate test for between-population selection
#'
#' Tests for the presence of directional and stabilizing selection on
#' multivariate traits for a multi-population dataset.
#'
#' @param G Genetic additive among-trait covariance matrix for k number of
#' traits. Can be either a single kxk matrix, a list of multiple Gs or
#' three-dimensional array (kxkxi), with i being the number of items.
#' @param means matrix containing the means/fixed effects of all k characters
#'     for each population (n). Can be either a single nxk matrix, a list of
#'     multiple nxk matrices or a three-dimensional array (nxkxi). Must match
#'     the dimensionality of G.
#' @param theta The coancestry coefficient among populations. Can be either a
#' single nxn matrix, a list of nxn matrices of a three-dimensional array
#' (nxnxi).
#' @param anc matrix containing the ancestral state of all k characters. Can
#' be either a single nxk matrix, a list of multiple nxk matrices or a
#' three-dimensional array (nxkxi). Must match the dimensionality of G. If not
#' provided, the ancestral state is assumed to be 0 for all traits.
#' @param Verbose logical. If TRUE, calculate the mahalanobis distance for each
#' population based on the expected divergence.
#' @param parallel Should be parallelized? Default is FALSE.
#' @return
#' If nsim=0 returns a single value for the parametric S-test. If nsim>0 returns
#' a non-parametric S-test value for each population.
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom plyr alply laply
#' @importFrom foreach foreach
#' @references
#' @author Fabio Andrade Machado
#'
driftsel<-function(G, means, theta, anc=NULL, verbose=TRUE, parallel=FALSE){
  if(is.null(anc)) {
    if(any(class(means)=="list")) {
      anc<-means
      for(i in seq_along(anc)) anc[[i]][]<-0
    } else {
      anc<-matrix(0, nrow(means),ncol(means))
    }
  }
  pars<-list(G=G, means=means, theta=theta, anc=anc)
  sizes<-foreach(i=seq_along(pars),.combine = "c") %do% {
    x<-pars[[i]]
    xname<-names(pars)[i]
    if(any(class(x)=="list")) {
      size<-length(x)
    } else {
      if(any(class(x)=="array")) {
        if(length(dim(x))==3) {
          size<-dim(x)[3]
        } else {
          size<-1
        }
      } else {
        stop(paste(xname, 'must be either a matrix, a cubic array or a list'))
      }
    }
    size
  }
  if(!all(sizes==max(sizes)|sizes==1)) stop('Number of iterations must be equal among G, means, theta and anc, or equal to 1')
  iters<-max(sizes)
  pars<-lapply(pars, function(x){
    if (any(class(x)!="list")){
      if (length(dim(x))==3) {
        x <- alply(x,3,identity)
      } else {
        l<-vector("list", iters)
        for(i in 1:iters) l[[i]]<-x
        l
      }
    } else x
  })

  D2 <- vector("numeric",iters)
  n<-nrow(pars$means[[1]])
  if(verbose) {
    D2_pop<-cdf_pop<-matrix(0,iters,n,dimnames = list(NULL,rownames(pars$means[[1]])))
    }

  for(i in seq_len(iters)){
    k<-nrow(pars$G[[i]])
    means <- pars$means[[i]]
    v_means <- c(means)
    Sigma <- 2*pars$G[[i]] %x% pars$theta[[i]]
    mu <- c(pars$anc[[i]])
    invSigma <- solve(Sigma)
    v_means <- mu - v_means
    D2[i]<-sum(v_means %*% invSigma * v_means)

    if(verbose){
      D2_pop[i,]<-laply(1:n, function(j){
        v_mean<-pars$anc[[i]][j,]-means[j,]
        sigma <- 2* pars$G[[i]] %x% pars$theta[[i]][j,j]
        invsigma <- solve(sigma)
        sum(v_mean %*% invsigma * v_mean)
      })
    }
  }
  cdf = pchisq(D2, df = k*n)
  cdf_pop[] = pchisq(D2_pop, df = k)

  if(verbose) {
    return(list(cdf=cdf, cdf_pop=cdf_pop, D2=D2, D2_pop=D2_pop))
    } else return(cdf)
}
