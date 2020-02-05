#' @name psc
#' @description Phylogenetic serial contrasts
#' @param x a numeric vector.
#' @param phy an object of class "phylo"
#' @return
#' lala
#' @importFrom ape pic ace reorder.phylo
#' @export
psc<-function(x, phy){
  allTraits<-c(x,ace(x,phy,method = "pic")$ace)
  rescaled.phy<-
    pic(rep(1,times=length(phy$tip.label)),phy,rescaled.tree = T)$rescaled.tree
  rescaled.phy<-
    reorder(rescaled.phy,order = "cladewise")
  seqc<-allTraits[phy$edge[,1]]-allTraits[phy$edge[,2]]
  seqc<-seqc/sqrt(rescaled.phy$edge.length)
  return(seqc)
}
