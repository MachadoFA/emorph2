#' sim_multiphylo simulates multivariate evolution.
#'
#' sim_multiphylo takes in a phylogeny and an aditive genetic
#' variance-covariance matrix and simulates trait evolution under different
#' selective regimes.
#'
#' @param phy a phylogenetic tree. Must be of class 'phylo'.
#' @param G matrix kxk for k number of characters.
#' @param n.s a vector indicating the sample sizes for eaxh terminal.
#' @param selection define the type oif selection "random" or numeric between 1
#'     and k.
#' @param efsize strenght of selection relative to drif. Must be >= 0.
#' @param gen_time generation time in the same time unity as the phylogenetic
#'     tree (usually MY).
#' @param Nef effective population size for initial population. Must be >= 1.
#' @param Nef_osc oscilation in effective population size (Nef) over time.
#'     Character "no", "norm", or "unif".
#' @param Nef_par set parameters for the oscilation of effective population
#'     size (Nef).
#' @param scale logic TRUE of FALSE. Should the matrices be scaled to the
#'     empirical character values (means)?
#' @param means matrix sxk containing the empirical averages of all k characters for each
#'     species (s). Default = NULL.
#' @param matrix logic TRUE of FALSE. Should the output be a set
#' of marices? See 'Value' for more details.
#'
#' @return
#' \describe{
#' The 'sim_multiphylo' function returns an object of class "list". This is a
#' list of items of class "matrix", all kxk.
#' If 'matrix' = TRUE, the function returns 5 variance covariance matrices:
#'   \item{R}{Evolutionary rate matrix. Describes the rate of evolution and
#'   co-evolution among characters.
#'   Diagonal contains traits' rate evolution according to Brownian-Motion,
#'   off diagonals represent traits co-evolution.}
#'   \item{W}{Within matrix. Pooled-within species' trait variance-covariance matrix. }
#'   \item{B}{Between matrix. Variance-vovariance between average species traits.}
#'   \item{G}{Genetic matrix. The original aditive variance-covariance G matrix imputed.}
#'   \item{A}{Selection matrix. The expected effect due to selection.}
#' If 'matrix' = FALSE, the function returns 3 matrices:
#'   \item{bdata}{Simulated trait values for species averages.}
#'   \item{wdata}{Simulated intraspecific error. }
#'   \item{G}{the original aditive variance-covariance G matrix imputed.}
#' }
#'

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
      svec<-matrix(rnorm(dim(G)[1]*selection), dim(G)[1], selection) %>% apply(., 2, evolqg::Normalize)
      C <- svec %*% t(svec) + diag(dim(G)[1]) * 0.0001
    }
    C <- C/tr(C)
    A <- (G %*% C %*% G) / tr(G %*% C %*% G)
    A <- efsize * A
    sel  <- phytools::sim.corrs(phy, A)
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
  drift<- phytools::sim.corrs(phy1, G)

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
    dplyr::group_by(.,sps) %>% dplyr::summarize_all(funs(mean))
  data.s<-data.s+de[,-1]

  if(matrix){
    pics<-apply(data.s, 2, function(x) pic(x, phy))
    out<-list(R=t(pics) %*% pics,
              W=var(ws),
              B=var(data.s),
              G=G,
              A=A)
  } else {
    out<-list(bdata=data.s,
              wdata=ws,
              G=G)
  }
  out
}


tr <- function(x) {sum(diag(x))}
