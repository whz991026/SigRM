#' Simulate data function
#'
#' @param n_Sites number of site
#' @param test_num number of sample in correct treatment
#' @param control_num number of sample in false treatment
#' @param min_expression min expression for abundance
#' @param max_expression max expression for abundance
#' @param per_me percentage of difference
#' @param prange1 parameter for simulate p1
#' @param p_s_1 parameter for step of p1
#' @param prange2 parameter for simulate p2
#' @param prange3 parameter for simulate p3
#' @param p_s_3 parameter for step of p3
#'
#' @return simulate read counts( meth read counts of control data,
#' meth read counts of test data, unmeth read counts of control data,
#' unmeth read counts of test data)
#'
#' @importFrom stats rbinom rnorm runif
#' @export
#'
#' @examples
#'  data <- simulateData (test_num=10,control_num=30)
#'
#'
simulateData<- function(n_Sites=2000,test_num=100,control_num=300,min_expression=1,max_expression=4,
                        per_me=1/2,prange1=c(0,0.2),
                        p_s_1=0.01, prange2=c(0.25,0.3),prange3=c(0,0.3),p_s_3=0.01){
  # Base expression
  q1 <- matrix(10^runif(n_Sites*control_num,min_expression,max_expression),
               n_Sites,control_num)
  q2 <- matrix(10^runif(n_Sites*test_num,min_expression,max_expression),
               n_Sites,test_num)


  # number of sites set to be modification site
  n_new <- round(n_Sites*(per_me))

  # set p value
  step <- ceiling((prange1[2]-prange1[1])/p_s_1)
  p_step <- round(n_new /step)

  if (ceiling((prange1[2]-prange1[1])/p_s_1)!=floor((prange1[2]-prange1[1])/p_s_1)){
    n_step <- floor((prange1[2]-prange1[1])/p_s_1)
  }else {
    n_step <- floor((prange1[2]-prange1[1])/p_s_1)-1
  }

  p_1 <-matrix(1,1,1)
  if (n_step !=0){
    for (i in 1 :n_step){
      p_1_ <- matrix(runif(p_step *control_num,p_s_1*(i-1)+prange1[1],p_s_1*i+prange1[1])
                     ,p_step,control_num)
      if (dim(p_1)[1]==1&dim(p_1)[2]==1){
        p_1 <- p_1_
      }else{
        p_1 <- rbind(p_1,p_1_)
      }
    }
    p_1_ <- matrix(runif((( n_new - p_step *n_step)*control_num),p_s_1*
                           n_step+prange1[1],prange1[2]),ncol=control_num)

    p_1 <- rbind(p_1,p_1_)

  }else{
    p_1 <- matrix(runif((( n_new - p_step *n_step)*control_num),p_s_1*
                          n_step+prange1[1],prange1[2]),ncol=control_num)
  }


  p_2 <- matrix(runif(round((n_Sites*(per_me))*test_num),prange2[1],prange2[2]),
                round((n_Sites*(per_me))),test_num)
  #a <- apply(p_1,1,mean)
  # index <- which(a<=0.02)
  #p_2[index,] <- matrix(runif(length(index)*test_num,4,5),
  #              length(index),test_num)*apply(p_1[index,],1,mean)

  n_new <- round(n_Sites*(1-per_me))
  step <- ceiling((prange3[2]-prange3[1])/p_s_3)
  p_step <- round(n_new /step)

  if (ceiling((prange3[2]-prange3[1])/p_s_3)!=floor((prange3[2]-prange3[1])/p_s_3)){
    n_step <- floor((prange3[2]-prange3[1])/p_s_3)
  }else {
    n_step <- floor((prange3[2]-prange3[1])/p_s_3)-1
  }

  p_3_1 <-matrix(1,1,1)
  p_3_2 <-matrix(1,1,1)
  if (n_step !=0){
    for (i in 1 :n_step){
      p_3_1_ <- matrix(runif(p_step *control_num,p_s_3*(i-1)+prange3[1],p_s_3*i+prange3[1])
                       ,p_step,control_num)
      p_3_2_ <- matrix(runif(p_step *test_num,p_s_3*(i-1)+prange3[1],p_s_3*i+prange3[1]),
                       p_step,test_num)
      if (dim(p_3_1)[1]==1&dim(p_3_1)[2]==1){
        p_3_1 <- p_3_1_
        p_3_2 <- p_3_2_
      }else{
        p_3_1 <- rbind(p_3_1,p_3_1_)
        p_3_2 <- rbind(p_3_2,p_3_2_)
      }
    }
    p_3_1_ <- matrix(runif((( n_new - p_step *n_step)*control_num),p_s_3*
                             n_step+prange3[1],prange3[2]),ncol=control_num)
    p_3_2_ <- matrix(runif(( n_new - p_step *n_step)*test_num,p_s_3*
                             n_step+prange3[1],prange3[2]),ncol=test_num)

    p_3_1 <- rbind(p_3_1,p_3_1_)
    p_3_2 <- rbind(p_3_2,p_3_2_)

  }else{
    p_3_1 <- matrix(runif((( n_new - p_step *n_step)*control_num),p_s_3*
                            n_step+prange3[1],prange3[2]),ncol=control_num)
    p_3_2 <- matrix(runif(( n_new - p_step *n_step)*test_num,p_s_3*
                            n_step+prange3[1],prange3[2]),ncol=test_num)
  }




  # p_1 is the control detection level
  p_1 <-rbind(p_1,p_3_1)
  # p_2 is the test detection level
  p_2 <-rbind(p_2,p_3_2)

  # simulate data
  meth1=matrix(1,n_Sites,control_num)
  meth2=matrix(1,n_Sites,test_num)
  unmeth1=matrix(1,n_Sites,control_num)
  unmeth2=matrix(1,n_Sites,test_num)


  for (j in 1: control_num){
    meth1[,j]=rbinom(n_Sites,round(q1)[,j],p_1[,j])
    unmeth1[,j]=rbinom(n_Sites,round(q1)[,j],(1-p_1)[,j])
  }

  for (j in 1: test_num){
    meth2[,j]=rbinom(n_Sites,round(q2)[,j],p_2[,j])
    unmeth2[,j]=rbinom(n_Sites,round(q2)[,j],(1-p_2)[,j])
  }


  res=list(4)
  res[[1]]=meth1 # control (num of detected)
  res[[2]]=meth2 # test (num of detected)
  res[[3]]=unmeth1 # control (num of undetected)
  res[[4]]=unmeth2 # test (num of undetected)

  return (res)
}
