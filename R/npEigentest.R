#' Non-parametric confidence intervals for eigen-decomposition-based tests of
#' matrix proportionality.
#'
#' @param G square matrix of covariance among for k number of traits.
#' @param phy a phylogenetic tree. Must be of class 'phylo'.
#' @param n.s a vector indicating the sample sizes for each terminal.
#' @param dim.ret number of dimensions that will be retained. Can be used to
#' remove small eigeinvalues.
#' @param parallel Should the analysis be parallelized? Default is FALSE. See
#' 'parallel vignette' for details.
#'
#' @return
#'  \describe{
#'  npEigentest returns a list containing summary of test results, simmulated
#'  values for each test, and observed results for empirical data.
#'  \itemize{
#'      \item{TestResults}{ summary table containing observed values and test results. TRUE = rejected drift.}
#'      \item{Empirical}{ Test results for empirical data.}
#'      \item{SimValues}{ simulated values for each test.}
#'      }
#'
#' @examples
#' data("Canidae")
#' test.out<-npEigentest(G=W, means, tree, n.s, sims = 100, dim.ret = 20)
#' test.out$SimValues
#' @export



npEigentest<-function(G,means,phy,n.s,sims=1000,dim.ret=NULL,parallel=FALSE){
  #Estimate the rate matrix from Independant contrasts. V/CV matrix of
  #evolutionary responses (DeltaZ per species) standardized by branch lenght
  #(divergence time).
  pics<-apply(means, 2, function(x) ape::pic(x, phy))
  R=t(pics) %*% pics

  # Empirical test
  obs<-pcTests(G,R,length(n.s))

  # simmulations
  sim<-adply(1:sims, 1, function(i){
    W <- mvtnorm::rmvnorm(sum(n.s),sigma=G) %>% var
    picsr <- mvtnorm::rmvnorm(dim(pics)[1],sigma=G)
    R <- t(picsr) %*% picsr/dim(pics)[1]
    pcTests(W,R,length(n.s))
  },.parallel = parallel)[,-1]

  quantiles <- sim %>% apply(.,2,function(c) stats::quantile(c, c(0.025,0.975)))


  # Quantiles for each test
  slt05<-obs$slt<quantiles[1,1] | obs$slt>quantiles[2,1] # if TRUE: reject drift
  SDrel05<-obs$SDrel<quantiles[1,2] | obs$SDrel>quantiles[2,2] # if TRUE: reject drift

  return(list("TestResults" = data.frame(obs,
                                  slt05,
                                  SDrel05),
              "Empirical" = obs,
              "SimValues" = sim
              ) )
}


pcTests<-function(G,R,n,dim.ret=NULL){
  eigenG<-eigen(G)
  if(is.null(dim.ret)) dim.ret=dim(R)[1]
  RGr<- t(eigenG$vectors) %*% R %*% eigenG$vectors
  df.s  <- data.frame(G=log(eigenG$value),
                      R=log(diag(RGr)))
  sltest<-lm(R~G, data=df.s[1:dim.ret,])
  sltest<-sltest$coefficients[2]

  CRGr<-stats::cov2cor(RGr)[1:dim.ret,1:dim.ret]
  dimnames(CRGr)<-list(1:dim(CRGr)[1],1:dim(CRGr)[1])
  evs<-eigen(CRGr)$values
  N<-length(evs)
  SDrel<-sqrt(sum((mean(evs)-evs)^2))/sqrt(N-1)

  if(n<dim.ret) {
    CRGr<-CRGr[1:n,1:n]
    evs<-eigen(CRGr)$values
    N<-length(evs)
    SDrel<-sqrt(sum((mean(evs)-evs)^2))/sqrt(N-1)
  }
  return(data.frame(slt=sltest[1],
                    SDrel))
}
