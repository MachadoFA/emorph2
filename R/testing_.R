sample_rg<-501:1000
afm$theta<-afm$theta[,,sample_rg]
afm$kappa<-afm$kappa[,,sample_rg]
afm$fst<-afm$fst[sample_rg]
afm$alpha<-afm$alpha[,sample_rg]

dims<-dim(Males_Fixed)[2]
full_Fixed<-array(NA,dim=c(12,dims*2,500),
                  dimnames = list(dimnames(Males_Fixed)[[1]],
                                  paste0(rep(c("f","m"),each=dims),
                                         dimnames(Males_Fixed)[[2]]),NULL))
full_Fixed[,1:dims,]<-Females_Fixed
full_Fixed[,1:dims+dims,]<-Males_Fixed
full_Gext<-full_G
iF<-0.25/2
scale<-1/(2*iF)
full_G<-full_G*scale


theta<-afm$theta[,,1]
fst<-afm$fst[1]
G<-full_G[,,1]
means<-full_Fixed[,,1]

gelli<-ellipse::ellipse(2*G*fst/(1-fst))

ggplot(data.frame(means),aes(fPC1,fPC2))+
  coord_fixed()+
  geom_polygon(aes(PC1.F,PC2.F),data.frame(gelli), alpha=0.3)+
  geom_point()

Sigma<-2*G%x%th

library(ggplot2)
library(mvtnorm)
library(plyr)
library(dplyr)

meansV<-rmvnorm(10000,sigma = Sigma)
x<-rowMeans(meansV)
dellis<-adply(meansV,1,function(x){
  x<-matrix(x,dim(means)[1],dim(means)[2])
  D <- t(x) %*% x / dim(means)[1]
  ellipse::ellipse(D)
})

meansV<-rmvnorm(10000,sigma = Sigma)
x<-colMeans(meansV)
x<-matrix(x,dim(means)[1],dim(means)[2])
D1 <- t(x) %*% x / dim(means)[1]

ggplot(data.frame(means),aes(fPC1,fPC2))+
  coord_fixed()+
  # geom_polygon(aes(x,y, group=X1),dellis[1:500,], alpha=0,color="black")+
  geom_polygon(aes(x,y),data.frame(ellipse::ellipse(D2)), alpha=0,color="black")+
  geom_point()



x<-means
cov<-G

foreach(i=1:500, .combine = "rbind") %do% {
  cov<-full_G[,,i]
  x<-full_Fixed[,,i]


  init<-Sys.time()
  mahalanobis(x,center=F, cov)
  end<-Sys.time()
  t1<-end-init

  init<-Sys.time()
  diag(x%*%solve(G)%*%t(x))
  end<-Sys.time()
  t2<-end-init

  init<-Sys.time()
  setNames(rowSums(x %*% solve(cov) * x), rownames(x))
  end<-Sys.time()
  t3<-end-init

  data.frame(t1,t2,t3)
} %>% reshape::melt() %>%
  ggplot(., aes(variable, as.numeric(value)))+
  geom_violin()+
  scale_y_log10()






