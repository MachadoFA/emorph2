library(testthat)

test_that("Test sim_multiphylo()", {
  tree<-readRDS(file="tree.RDS")
  G<-readRDS(file="G.RDS")
  ns<-readRDS(file="ns.RDS")

  expect_equal(matrix(0,2,2), sim_multiphylo(G = G,phy = tree,n.s = ns,efsize=0)$A)
})

