#' Fix and fuse replicates
#'
#'
#' @param rep1
#' @param rep1
#' @param pairdLM
#'
#' @importFrom Morpho fixLMmirror
#' @importFrom Morpho fixLMtps
#' @importFrom shapes procOPA
#' @author Fabio Andrade Machado

fixNfuseReps<-function(rep1, rep2, pairedLM){
  for(i in 1:dim(rep1)[3]){
    x<-abind(rep1[,,i],rep2[,,i]) %>%
      fixLMmirror(.,pairedLM = pairedLM)
    miss<-!is.na(x[,1,1])
    xx<-fixLMtps(x[miss,,])
    opa<-procOPA(xx[,,1],xx[,,2], scale=F)
    xx<-mshape(abind(opa$Ahat,opa$Bhat))
    x[miss,,]<-xx
    x<-x[,,1]
    if(i==1) coords_tmp<-x else coords_tmp<-abind(coords_tmp,x)
  }
  rownames(coords_tmp)<-rownames(rep1)
  coords_tmp
}
