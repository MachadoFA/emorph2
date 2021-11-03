#' Variance of the product a random variable by a scalar
#'
#' @param var Variance of a random variable
#' @param s multiplying scalar
#'
#' @author Fabio Andrade Machado
#' @export
varScalar<-function(var,s){
  var*s^2
}

#' Variance of the sum of two random variables
#'
#' @param V Covariance matrix of two traits
#'
#' @author Fabio Andrade Machado
#' @export
varSum<-function(V){
  var1<-V[1,1]
  var2<-V[2,2]
  cov <-V[1,2]
  var1+var2+2*cov
}

#' Variance of the product of two random variables
#'
#' @param V Covariance matrix of two traits
#' @param m Vector of trait means
#'
#' @author Fabio Andrade Machado
#' @export
varProd<-function(V,m){
  varx<-V[1,1]
  vary<-V[2,2]
  cov <-V[1,2]
  mx  <-m[1]
  my  <-m[2]

  mx^2*vary + my^2*varx + cov^2 + 2*mx*my*cov + varx*vary
}

#' Variance of the ratio of two random variables
#'
#' @param V Covariance matrix of two traits
#' @param m Vector of trait means
#'
#' @author Fabio Andrade Machado
#' @export
varRat<-function(V,m){
  varx<-V[1,1]
  vary<-V[2,2]
  cov <-V[1,2]
  mx  <-m[1]
  my  <-m[2]

  varx/my^2 + mx^2*vary/my^4 - 2*mx*cov/my^3
}

#' Expectation of the ratio of two random variables
#'
#' @param V Covariance matrix of two traits
#' @param m Vector of trait means
#'
#' @author Fabio Andrade Machado
#' @export
expRat<-function(V,m){
  varx<-V[1,1]
  vary<-V[2,2]
  cov <-V[1,2]
  mx  <-m[1]
  my  <-m[2]

  mx/my - cov/my^2 +vary*mx/my^3
}


#' Covariance of the product a random variable by a scalar
#'
#' @param cov Covariance between two random variables
#' @param s1 multiplying scalar for the first variable
#' @param s2 multiplying scalar for the second variable
#'
#' @author Fabio Andrade Machado
#' @export
covScalar<-function(cov,s1=1,s2=1){
  cov*s1*s2
}

#' Covariance of a random variable and the sum of two random variables
#'
#' @param V Covariance matrix of three traits, the second and third traits are the ones being summed.
#' @param covxy Covariance between two random variables
#' @param covxu Covariance between two random variables
#'
#' @author Fabio Andrade Machado
#' @export
covSum<-function(V, covxy=NULL, covxu=NULL){
  if(is.null(covxy)|is.null(covxu)){
    covxy<-V[1,2]
    covxu<-V[1,3]
  }
  covxy + covxu
}

#' Covariance of a random variable and the product of two random variables
#'
#' @param V Covariance matrix of three or four traits.
#' @param m Vector of expectations of the variables being multiplied
#'
#' @author Fabio Andrade Machado
#' @export
covProd<-function(V, m){
  if(all(dim(V)==3)){
    covyv<-V[2,3]
    covxv<-V[1,3]

    mx<-m[1]
    my<-m[2]

    pvar<-mx*covyv + my*covxv
  }
  if(all(dim(V)==4)){
    covyv<-V[2,4]
    covyu<-V[2,3]
    covxv<-V[1,4]
    covxu<-V[1,3]

    mx<-m[1]
    my<-m[2]
    mu<-m[3]
    mv<-m[4]

    pvar<-mx*mu*covyv + mx*mv*covyu + my*mu*covxv + my*mv*covxu + covxu*covyv + covxv+covyu
  }
  return(pvar)
}

