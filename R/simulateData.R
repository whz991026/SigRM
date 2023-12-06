#' Simulate data function
#'
#' @param n_Sites number of site
#' @param test_num number of sample in correct treatment
#' @param control_num number of sample in false treatment
#' @param min_expression min expression for abundance
#' @param max_expression max expression for abundance
#' @param expression_p expression parameter
#' @param per_me percentage of difference
#' @param prange1 parameter for simulate p1
#' @param prange2 parameter for simulate p2
#' @param prange3 parameter for simulate p3
#' @param shape2 the hyper-parameter of the beta distribution
#'
#' @return simulate read counts( meth read counts of control data,
#' meth read counts of test data, unmeth read counts of control data,
#' unmeth read counts of test data)
#'
#' @importFrom stats rbinom rnorm runif rbeta
#' @export
#'
#' @examples
#'  data <- simulateData (test_num=10,control_num=30)
#'
#'
simulateData <- function(n_Sites=2000,test_num=100,control_num=300,min_expression=0,max_expression=4,
         expression_p=10,per_me=1/2,prange1=c(0,0.4),prange2=c(0.6,0.8),
         prange3=c(0,0.4),shape2=5){
  # Base expression
  q1 <- matrix(expression_p^runif(n_Sites*control_num,min_expression,max_expression),
               n_Sites,control_num)
  q2 <- matrix(expression_p^runif(n_Sites*test_num,min_expression,max_expression),
               n_Sites,test_num)
  
  
  # number of sites set to be modification site
  n_new <- round(n_Sites*(per_me))
  
  p_1 <- as.list(runif(n_new,prange1[1],prange1[2]))
  p_1 <- t((as.data.frame(lapply(p_1, function(x) {
    shape1 <- shape2*x/(1-x)
    return(rbeta(control_num,shape1,shape2))
  }))))
  rownames(p_1) <-c()
  
  p_2 <- matrix(runif(n_new*test_num,prange2[1],prange2[2]),ncol=test_num)
  
  
  p_3 <- as.list(runif((n_Sites-n_new),prange3[1],prange3[2]))
  p_3 <- t(as.data.frame(lapply(p_3, function(x) {
    shape1 <- shape2*x/(1-x)
    return(rbeta((control_num+test_num),shape1,shape2))
  })))
  rownames(p_3) <-c()
  # simulate data
  
  p_1 <- rbind(p_1,p_3[,1:control_num])
  p_2 <- rbind(p_2,as.matrix(p_3[,-(1:control_num)]))
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(q1)[1],dim(q1)[2]),q1,p_1)))
  list_function <- function(x){
    meth_control=rbinom(x[1],round(x[2:(x[1]+1)]),x[(2+x[1]):(2*x[1]+1)])
  }
  meth_control <- as.data.frame(lapply(list_create,list_function))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(q1)[1],dim(q1)[2]),q1,p_1)))
  list_function <- function(x){
    unmeth_control=rbinom(x[1],round(x[2:(x[1]+1)]),(1-x[(2+x[1]):(2*x[1]+1)]))
  }
  unmeth_control <- as.data.frame(lapply(list_create,list_function))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(q2)[1],dim(q2)[2]),q2,p_2)))
  list_function <- function(x){
    meth_test=rbinom(x[1],round(x[2:(x[1]+1)]),x[(2+x[1]):(2*x[1]+1)])
  }
  meth_test <- as.data.frame(lapply(list_create,list_function))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(q2)[1],dim(q2)[2]),q2,p_2)))
  list_function <- function(x){
    unmeth_test=rbinom(x[1],round(x[2:(x[1]+1)]),(1-x[(2+x[1]):(2*x[1]+1)]))
  }
  unmeth_test <- as.data.frame(lapply(list_create,list_function))
  
  res=list(4)
  res[[1]]=meth_control # control (num of detected)
  res[[2]]=meth_test # test (num of detected)
  res[[3]]=unmeth_control # control (num of undetected)
  res[[4]]=unmeth_test # test (num of undetected)
  
  return (res)
}