#' Shifting nodes
#'
#' rescaling phy
#'
#' @param phy an object of class "phylo"
#' @param criteria icriteri
#' @details S
#' @return A
#' @importFrom ape pic ace reorder.phylo
#' @export
#' @author Fabio Andrade Machado
node_shift<-function(phy,criteria,shift=NULL){
  phy<-reorder(phy,order = "postorder")
  if(is.null(shift)) shift<-criteria

  for(i in 1:nrow(phy$edge)){
    if(phy$edge.length[i]<criteria){
      node2shift <- phy$edge[i,1]
      phy$edge.length[phy$edge[,2]==node2shift]<-phy$edge.length[phy$edge[,2]==node2shift]-shift
      phy$edge.length[phy$edge[,1]==node2shift]<-phy$edge.length[phy$edge[,1]==node2shift]+shift
    }
  }

  phy<-reorder(phy)

  for(i in 1:nrow(phy$edge)){
    if(phy$edge.length[i]<criteria){
      node2shift <- phy$edge[i,1]
      phy$edge.length[phy$edge[,2]==node2shift]<-phy$edge.length[phy$edge[,2]==node2shift]-shift
      phy$edge.length[phy$edge[,1]==node2shift]<-phy$edge.length[phy$edge[,1]==node2shift]+shift
    }
  }


  phy
}
