#' npEigentest produces non-parametric confidence intervals spectral
#' decomposition-based tests of matrix proportionality.
#'
#' npEigentest yadayaadayda.
#'
#' @param G matrix kxk for k number of traits.
#' @param phy a phylogenetic tree. Must be of class 'phylo'.
#' @param n.s a vector indicating the sample sizes for each terminal.

npEigentest<-function(G,means,phy,n.s,sims=1000,dim.ret=NULL,parallel=FALSE){
  pics<-apply(means, 2, function(x) pic(x, phy))
  R=t(pics) %*% pics

  obs<-pcTests(G,R,length(n.s))

  sim<-adply(1:sims, 1, function(i){
    W <- rmvnorm(sum(n.s),sigma=G) %>% var
    picsr <- rmvnorm(dim(pics)[1],sigma=G)
    R <- t(picsr) %*% picsr/dim(pics)[1]
    pcTests(W,R,length(n.s))
  },.parallel = parallel)[,-1] %>% apply(.,2,function(c) quantile(c, c(0.025,0.975)))

  slt05<-obs$slt<sim[1,1] | obs$slt>sim[2,1]
  SDrel05<-obs$SDrel<sim[1,2] | obs$SDrel>sim[2,2]

  data.frame(obs,slt05,SDrel05)
}


pcTests<-function(G,R,n,dim.ret=NULL){
  eigenG<-eigen(G)
  if(is.null(dim.ret)) dim.ret=dim(R)[1]
  RGr<- t(eigenG$vectors) %*% R %*% eigenG$vectors
  df.s  <- data.frame(G=log(eigenG$value),
                      R=log(diag(RGr)))
  sltest<-lm(R~G, data=df.s[1:dim.ret,])
  sltest<-sltest$coefficients[2]

  CRGr<-cov2cor(RGr)[1:dim.ret,1:dim.ret]
  dimnames(CRGr)<-list(1:dim(CRGr)[1],1:dim(CRGr)[1])
  evs<-eigen(CRGr)$values
  N<-length(evs)
  SDrel<-sqrt(sum((mean(evs)-evs)^2))/sqrt(N-1)

  if(n<dim.ret) {
    CRGr<-CRGr[1:n,1:n]
    evs<-eigen(CBGr)$values
    N<-length(evs)
    SDrel<-sqrt(sum((mean(evs)-evs)^2))/sqrt(N-1)
  }
  data.frame(slt=sltest[1],
             SDrel)
}
