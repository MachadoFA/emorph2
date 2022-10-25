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
#' @param sims Numeric. Number of simulations used for the non-parametric tests.
#' If NULL, the parametric test is performed.
#' @param pop.test logical. If TRUE, calculate the mahalanobis distance for each
#' population based on the expected divergence.
#' simulations.
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
driftsel<-function(G, means, theta, anc=NULL, sims=0, parallel=FALSE){
  if(is.null(anc)) {
    if(any(class(means)=="list")) {
      anc<-matrix(0, dim(means[[1]])[1],dim(means[[1]])[2])
    } else {
      anc<-matrix(0, dim(means)[1],dim(means)[2])
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
  if (parallel) `%do_%`<-`%dopar%` else `%do_%`<-`%do%`
  pars<-lapply(pars, function(x){
    if (any(class(x)!="list")){
      if (length(dim(x))==3) {
        x <- alply(x,3,identity)
      } else {
        l<-vector("list", iters)
        for(i in 1:iters) l[[i]]<-x
        l
      }
    }
  })
  k<-nrow(pars$G[[1]])
  n<-nrow(pars$means[[1]])
  D2 <- vector("numeric",iters)
  if(sims>1) {
    cdf_pop<- matrix(0,iters,n,dimnames = list(NULL,rownames(pars$means[[1]])))
  }

  for(i in seq_len(iters)){
    means <- pars$means[[i]]
    v_means <- c(means)
    Sigma <- 2*pars$G[[i]] %x% pars$theta[[i]]
    mu <- c(pars$anc[[i]])
    invSigma <- solve(Sigma)
    v_means <- mu - v_means
    D2[i]<-rowSums(v_means %*% invSigma * v_means)

    if(sims>0){
      D2_pop<-laply(1:n, function(j){
        means[-j,]<-pars$anc[[i]][-j,]
        v_means <- mu - c(means)
        rowSums(v_means %*% invSigma * v_means)
      })
      x<-rmvnorm(sims,sigma = Sigma)
      D2_pop_sims<-foreach(j=seq_len(sims),.combine = "rbind") %do_% {
        means_sim<-matrix(x[j,],ncol = k,byrow = TRUE)
        laply(1:n, function(k){
          means_sim[-k,]<-pars$anc[[i]][-k,]
          v_means <- mu - c(means_sim)
          rowSums(v_means %*% invSigma * v_means)
        })
      }
      cdf_pop[i,]<-laply(1:n, function(j){
        mean(D2_pop[j]>D2_pop_sims[,j])
      })
    }
  }
  cdf = pchisq(D2, df = k*n)
  if(sims>0) return(list(cdf,cdf_pop)) else return(cdf)
}