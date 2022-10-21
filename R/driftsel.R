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
#' @param parallel Should be parallelized? Default is FALSE. See 'parallel vignette'
#' for details.
#' @return
#' If boot=FALSE and verbose=FALSE, returns a single value for the parametric
#' S-test.
#' If boot=FALSE and verbose=TRUE, returns a distribution of posterior
#' probabilities for the S-test.
#' If boot=TRUE and verbose=FALSE, returns the empirical and null distribution
#' of squared mahalanobis distances (D2).
#' If boot=TRUE and verbose=TRUE, in addition to the previous result, returns
#' the following
#'  \describe{
#'      \item{pop}{A list containing the empirical and bootstrapped mahalanobis
#'      distance for each population}
#'      \item{null}{The distribution of the expected among-population
#'      covariances according to genetic drift.}
#'      }
#' @export
#' @importFrom mvtnorm rmvnorm
#' @importFrom plyr alply laply
#' @importFrom foreach foreach
#' @references
#' @author Fabio Andrade Machado
#'
driftsel<-function(G, means, theta, anc=NULL, sims=1, pop.test=FALSE){

  if(is.null(anc)) if(any(class(means)=="list")) {
    anc<-matrix(0, dim(mean[[1]])[1],dim(mean[[1]])[2])
  } else {
    anc<-matrix(0, dim(means)[1],dim(means)[2])
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

  k<-nrow(pars$G[[1]])
  n<-nrow(pars$means[[1]])
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
    }
  })

  D2 <- vector("numeric",iters)
  if(pop.test & sims>1) {
    D2_pop <- matrix(0,iters,n,dimnames = list(NULL,colnames(pars$means[[i]])))
    D2_pop.tmp <- D2_pop
    D2_sims<-vector("list",iters)
    D2_pop_sims<-vector("list",iters)
  }

  for(i in seq_len(iters)){
    means <- pars$means[[i]]
    v_means <- c(means)
    Sigma <- 2*pars$G[[i]] %x% pars$theta[[i]]
    mu <- c(pars$anc[[i]])
    invSigma <- solve(Sigma)
    v_means <- mu - v_means
    D2[i]<-rowSums(v_means %*% invSigma * v_means)
    # t(mu - v_means) %*% invSigma %*% (mu - v_means)

    if(pop.test & sims>1){
      D2_pop[i,]<-laply(1:n, function(j){
        means[-j,]<-0
        v_means <- mu - c(means)
        rowSums(v_means %*% invSigma * v_means)
      })
    }

    if(sims>1){
      x<-rmvnorm(sims,sigma = Sigma)
      D2_sims[[i]]<-rowSums(x %*% invSigma * x)

      if(pop.test){
        D2_pop.tmp[i,]<-laply(1:n, function(j){
          means[-j,]<-0
          v_means <- mu - c(means)
          rowSums(v_means %*% invSigma * v_means)
        })
        D2_pop_sims[[i]]<-D2_pop.tmp
      }
    }
  }

  cdf = pchisq(D2, df = k*n)

  if (parallel) `%do_%`<-`%dopar%` else `%do_%`<-`%do%`

  # if(sims>1){
  #   foreach (i=seq_len(iters),.combine = "rbind") %do% {
  #     Sigma <- 2*pars$G[[i]] %x% pars$theta[[i]]
  #     # mu <- c(pars$anc[[i]])
  #     invSigma <- solve(Sigma)
  #     # v_means <- mu - v_means
  #     x<-rmvnorm(sims,sigma = Sigma)
  #     globalD2<-rowSums(x %*% invSigma * x)
  #   }
  # }



  #   D2sim <- foreach (i=seq_len(iters),.combine = "rbind") %do% {
  #     means <- pars$means[[i]]
  #     Sigma <- 2*pars$G[[i]] %x% pars$theta[[i]]
  #     # mu <- c(pars$anc[[i]])
  #     invSigma <- solve(Sigma)
  #     # v_means <- mu - v_means
  #     x<-rmvnorm(sims,sigma = Sigma)
  #     globalD2<-rowSums(x %*% invSigma * x)
  #
  #     x<-array(x,dim = c(n,k,sims))
  #     if(pop.test){
  #       popD2<-foreach(j=1:sims,.combine = "rbind") %do_% {
  #         if(pop.test){
  #           ED <- t(x[,,i]) %*% x[,,i] / n
  #           setNames(rowSums(means %*% solve(ED) * means), rownames(means))
  #         }
  #       }
  #       data.frame(global=globalD2,popD2)
  #     }
  #   }
  #   if(pop.test) {
  #     return(list(D2=D2, postprob=cdf, globalsims=D2sim[,1], pop.test=D2sim[,-1]))
  #     } else return(list(D2=D2, postprob=cdf, globalsims=D2sim[,1]))
  # } else return(list(D2=D2, postprob=cdf))
}
