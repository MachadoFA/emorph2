#' Non-parametric tests of matrix proportionality for comparative data
#'
#' Evaluates the proportionality between the intraspecific patterns of
#' covariation and interespecific rates of evolution through eigen-decompositon.
#'
#' @param G matrix kxk for k number of traits.
#' @param means matrix sxk containing the empirical means of all k characters
#'     for each species (s). Default = NULL.
#' @param phy a phylogenetic tree. Must be of class 'phylo'.
#' @param n.s a vector indicating the sample sizes for each terminal.
#' @param sims numeric. total number of summulations.
#' @param dim.ret choose the number of dimentions that should be retained in the
#'     analysis.
#' @param parallel Should parallelize? Default is FALSE. See 'parallel vignette'
#'     for details
#' @return
#'  \describe{
#'  npEigentest returns a list containing summary of test results, simulated
#'  values for each test, and observed results for empirical data.
#'  \itemize{
#'      \item{Empirical}{ Empirical statistic values.}
#'      \item{SimValues}{ Simulated statistic values.}
#'      \item{eigenvalues}{ Relationship between empirical within-species
#'      eigenvalues and between species variances.}
#'      }}
#' @examples
#' \dontrun{data("Canidae")}
#' \dontrun{test.out<-npEigentest(G=W, means, tree, n.s, sims = 100, dim.ret = 20)
#' test.out$SimValues}
#' @export
#' @importFrom ape pic
#' @importFrom mvtnorm rmvnorm
#' @importFrom plyr adply
#' @importFrom dplyr %>%

eigenRegr<-function(G,means,phy,n.s,sims=1000,dim.ret=NULL,parallel=FALSE){
  #Estimate the rate matrix from Independant contrasts. V/CV matrix of
  #evolutionary responses (DeltaZ per species) standardized by branch lenght
  #(divergence time).
  pics<-apply(means, 2, function(x) pic(x, phy))
  R=t(pics) %*% pics

  # Empirical test
  obs<-pcTests(G,R,length(n.s), dim.ret=dim.ret)
  # Obtaining empirical distribution of eigenvalues
  eigenG<-eigen(G)
  if(is.null(dim.ret)) dim.ret=dim(R)[1]
  RGr<- t(eigenG$vectors) %*% R %*% eigenG$vectors
  df.s  <- data.frame(G=log(eigenG$value)-mean(log(eigenG$value)),
                      R=log(diag(RGr))-mean(log(diag(RGr))))

  # simulations
  sim<- plyr::adply(1:sims, 1, function(i){
    W <- rmvnorm(sum(n.s),sigma=G) %>% var
    picsr <- rmvnorm(dim(pics)[1],sigma=G)
    R <- t(picsr) %*% picsr/dim(pics)[1]
    pcTests(W,R,length(n.s),dim.ret=dim.ret)
  },.parallel = parallel)[,-1]

  quantiles <- sim %>% apply(.,2,function(c) stats::quantile(c, c(0.025,0.975)))

  # Quantiles for each test
  # slt05<-obs$slt<quantiles[1,1] | obs$slt>quantiles[2,1] # if TRUE: reject drift
  # SDrel05<-obs$SDrel<quantiles[1,2] | obs$SDrel>quantiles[2,2] # if TRUE: reject drift

  return(list("Empirical" = obs,
              "SimValues" = sim,
              "eigenvalues" = df.s))
  }


pcTests<-function(G,R,n,dim.ret=NULL){
  eigenG<-eigen(G)
  if(is.null(dim.ret)) dim.ret=dim(R)[1]
  RGr<- t(eigenG$vectors) %*% R %*% eigenG$vectors
  df.s  <- data.frame(G=log(eigenG$value)-mean(log(eigenG$value)),
                      R=log(diag(RGr))-mean(log(diag(RGr))))
  sltest<-lm(R~G, data=df.s[1:dim.ret,])
  sltest<-sltest$coefficients

  # CRGr<-stats::cov2cor(RGr)[1:dim.ret,1:dim.ret]
  # dimnames(CRGr)<-list(1:dim(CRGr)[1],1:dim(CRGr)[1])
  # evs<-eigen(CRGr)$values
  # N<-length(evs)
  # SDrel<-sqrt(sum((mean(evs)-evs)^2))/sqrt(N-1)
  #
  # if(n<dim.ret) {
  #   CRGr<-CRGr[1:n,1:n]
  #   evs<-eigen(CRGr)$values
  #   N<-length(evs)
  #   SDrel<-sqrt(sum((mean(evs)-evs)^2))/sqrt(N-1)
  # }
  return(data.frame(intercept=sltest[1],
                    slope=sltest[2] #,SDrel
                    ))
}
