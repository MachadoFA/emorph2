########################################
#simulating phylogenetic data###########
########################################
require(phytools)
require(plyr); require(dplyr)
#require(psych)
#require(doParallel)
require(mvtnorm)

tr<-function(x) {sum(diag(x))}

#registerDoParallel(cores = round(detectCores()/2))

# G- Pooled within covariance matrix
# means- species means
# phy- phylogeny
# n.s- vector of sample sizes
# inter- number of interaction
# parallel- run function on parallel (recomended)

# a<-"266"
# G<- AllCovs$vcv[[a]]
# phy <- extract.clade(filogenia, a)
# means <- allmedias[phy$tip.label,]
# n.s<- unlist(AllCovs$df[rownames(means)])+1
# selection="random"
# efsize=0
# gen_time=1e-06
# Nef=10000
# Nef_osc="no"

sim_multiphylo <-function(G,
                          phy,
                          n.s,
                          selection="random",
                          efsize=0,
                          gen_time=1e-06,
                          Nef=1000,
                          Nef_osc="no",
                          matrix=TRUE,
                          Nef_par=NULL,
                          scale=FALSE,
                          means=NULL){

  phy$edge.length <- phy$edge.length / gen_time
  n<-length(phy$tip.label)

  G <- G / tr(G)

  if(efsize==0) {sel <- 0; A<-diag(rep(0, nrow(G)))} else {
    if(selection=="random") C <- diag(dim(G)[1])
    if(is.numeric(selection)) {
      svec<-NULL
      svec<-matrix(rnorm(dim(G)[1]*selection), dim(G)[1], selection) %>% apply(., 2, Normalize)
      C <- svec %*% t(svec) + diag(dim(G)[1]) * 0.0001
    }
    C <- C/tr(C)
    A <- (G %*% C %*% G) / tr(G %*% C %*% G)
    A <- efsize * A
    sel  <- sim.corrs(phy, A)
  }

  if(Nef_osc=="no")   Nef_osc <- Nef else{
    if(Nef_osc=="norm") Nef_osc <- Nef * exp(rnorm(length(phy$edge.length), sd=Nef_par)) else{
      if(Nef_osc=="unif") Nef_osc <- runif(length(phy$edge.length), min = Nef_par[1], max = Nef_par[2])
    }
  }

  phy1 <- phy
  phy1$edge.length <- phy$edge.length / Nef_osc
  phy$edge.length  <- phy$edge.length / Nef
  # set.seed(seed)
  drift<- sim.corrs(phy1, G)

  # browser()
  data.s <- (drift + sel)
  if(scale) {
    pics<-apply(data.s, 2, function(x) pic(x, phy))
    rate.s <- tr(t(pics) %*% pics)
    s      <- sqrt(rate/rate.s)
    data.s <- data.s * s
  }
  # set.seed(seed)
  ws<-rmvnorm(sum(n.s),sigma=G)
  sps<-rep(rownames(data.s),times=n.s)
  de<-data.frame(sps=factor(sps,unique(sps)),ws) %>%
    group_by(.,sps) %>% summarize_all(funs(mean))
  data.s<-data.s+de[,-1]

  if(matrix){
    pics<-apply(data.s, 2, function(x) pic(x, phy))
    out<-list(R=t(pics) %*% pics,
              W=var(ws),
              B=var(data.s),
              G=G,
              A=A)
  } else {
    out<-list(bdata=data.s,wdata=ws,G=G)
  }
  out
}
