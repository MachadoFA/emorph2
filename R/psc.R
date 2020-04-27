#' Phylogenetic serial contrasts
#'
#' Calculate McPeek's phylogenetic contrasts along all branches of a phylogeny
#'
#' @param x a numeric vector.
#' @param phy an object of class "phylo"
#' @param scale if contrasts should be scaled by expected variance. Default is
#' TRUE
#' @details Sequential contrasts are calculated along the branches of the
#' phylogeny, instead of being calculated at each node, like standard
#' Phylogenetic Independent Contrasts (Felsenstein, 1985). Like PIC, these
#' contrasts can be scaled by expected variance and are a localized measure of
#' evolutionary rates. These contrasts where first proposed by McPeek (1995) to
#' find branches on a phylogeny where rates of evolution where higher.
#' @return A vector of phylogenetically sequential contrasts. Vector contains
#' one entry for each branch of the tree on the order that they appear on the
#' phylo object (cladewise order. See \code{\link[ape]{reorder.phylo}}).
#' @importFrom ape pic ace reorder.phylo
#' @export
#' @references McPeek, M. 1995. "Testing Hypotheses About Evolutionary Change on
#' Single Branches of a Phylogeny Using Evolutionary Contrasts". The American
#' Naturalist 145(5):686–703.
#' @references Felsenstein, J. 1985. Phylogenies and the comparative method.
#' The American Naturalist 125(1):1–15.
#' @author Fabio Andrade Machado
#' @seealso \code{\link[ape]{pic}}
psc<-function(x, phy, scale=TRUE){
  allTraits<-c(x,ace(x,phy,method = "pic")$ace)
  seqc<-allTraits[phy$edge[,1]]-allTraits[phy$edge[,2]]
  if(scale){
    rescaled.phy<-
      pic(x,phy,rescaled.tree = T)$rescaled.tree
    rescaled.phy<-
      reorder(rescaled.phy,order = "cladewise")
    seqc<-seqc/sqrt(rescaled.phy$edge.length)
  }
  return(seqc)
}

