#' This is the document that contains the sub function for the main function



#' Calculate the size factors.
#'
#' @param data read counts
#'
#' @return size factor
#' @importFrom stats median
sizeFactor <- function(data) {
  # make the elements no smaller than 0
  data <- as.matrix(data)
  data <- pmax(data,1e-5)

  # first log
  log_data <- log(data)
  log_mean <- rowMeans(log_data)
  log_s <- log_data-log_mean

  # then exp
  s_size <- exp(apply(log_s,2,median))
  return(s_size)
}




#' estimate P function.
#'
#' @param meth detection read counts data
#' @param unmeth total read counts minus detection read counts data
#' @param size size factor
#'
#' @return an estimated P vector
estimateP <- function(meth, unmeth, size) {
  p <- rowSums(t(t(meth)/size))/rowSums(t(t(unmeth+meth)/size))
  p[is.na(p)] <- 0.5
  p[is.infinite(p)] <- 0.5
  return(p)
}

#' estimate Q function.
#'
#' @param total_control total read counts for control data
#' @param total_test total read counts for test data
#' @param s_control size factor for control data
#' @param s_test size factor for test data
#'
#' @return an estimated q vector
estimateQ <- function(total_control,total_test,s_control,s_test){
  temp <- t(t(cbind(total_control,total_test))/c(s_control,s_test))
  q <- rowMeans(temp)
  return(q)
}

#' estimate E function.
#'
#' @param meth detection read counts data
#' @param unmeth total read counts minus detection read counts data
#' @param size size factor
#' @param q estimate the abundance of feature base on all data/ false treatment data
#'
#' @return an estimated e vector
estimateE <- function(meth,unmeth,size,q){

  e <- rowMeans(t(t(meth+unmeth)/size))/q
  e[is.na(e)] <- 0
  return(e)
}

#' The function to get estimated W
#'
#' @param data the data of read count
#' @param size size factor
#' @param e e factor
#'
#' @return a estimated W
calculateW <- function(data,size,e){

  q <-  t(t(data)/size)/e
  q <- q-rowMeans(q)
  w <- rowSums(q^2)/(length(size)-1)
  w <- as.matrix(w)
  w <- pmax(w,1e-8)
  return(w)
}


#' The function to get estimated Z
#'
#' @param q estimated abundance vector
#' @param p the estimated p
#' @param size the estimated size factor
#' @param e the estimated e vector
#'
#' @return estimated Z
calculateZ <- function(q,p,size,e){
  temp <- p*q/length(size)

  z <- rowSums(1/e%*%t(size))*temp


  return(z)
}


#' The function to local fit the W
#'
#' @param p estimated p level
#' @param q estimated q abundance
#' @param w estimated W
#'
#' @return a locfit list
#' @import locfit
locfitW <- function(p,q,w) {
  l <- log(q+1)
  intercept<-c(rep(1,length(q)))
  data<-data.frame(cbind(p,l,intercept,w))
  fit=locfit(w~lp(p,l,intercept),data=data,family="gamma")
  return(fit)
}


#' The function to get final estimated Z
#'
#' @param p estimated p level
#' @param q estimated q abundance
#' @param fit the fit list in the locfitw function
#'
#' @return the final estimated W
#' @importFrom stats predict
fittedW <- function(p,q,fit){
  l <- log(q+1)
  intercept<-c(rep(1,length(q)))
  data=data.frame(cbind(p,l,intercept))
  w_fit <- predict(fit,data)
  return(w_fit)
}



#' The function to get the fit of W
#'
#' @param meth detection read counts data
#' @param unmeth total read counts minus detection read counts data
#' @param p estimated p vector
#' @param q estimated q vector
#' @param e estimated e vector
#' @param size size factor
#'
#' @return a list of the fit variable of meth data and unmeth data
baseFit <- function(meth,unmeth,p,q,e,size){
  # get estimated w
  w_meth <-calculateW(meth,size,e)
  w_unmeth <-calculateW(unmeth,size,e)

  # locfit
  fit_meth <- locfitW(p,q,w_meth)
  fit_unmeth <- locfitW(p,q,w_unmeth)

  res <- list(fit_meth,fit_unmeth)
  return(res)
}

#' The function to calculate all estimate variable.
#'
#' @param meth_control detect read counts from false treatment
#' @param meth_test detect read counts from correct treatment
#' @param unmeth_control total read counts minus detection read counts data from false treatment
#' @param unmeth_test total read counts minus detection read counts data from correct treatment
#' @param s_control size factor for control data
#' @param s_test size factor for test data
#'
#'
#'
#' @return a list
getEstimate <- function(meth_control,meth_test,unmeth_control,unmeth_test,s_control,s_test){

  # estimate probability of methylation under a condition
  p_control <- estimateP(meth_control, unmeth_control, s_control)
  if (length(dim(meth_test)[2])==0){

    p_test <- estimateP(meth_test, unmeth_test, s_test)
    p0 <- estimateP(cbind(meth_control,meth_test), cbind(unmeth_control,
                                                         unmeth_test), c(s_control,s_test))

  }else{p_test <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
  p0 <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
  for (i in 1:dim(meth_test)[2]){
    p_test[,i] <- estimateP(meth_test[,i], unmeth_test[,i], s_test[i])
    p0[,i] <- estimateP(cbind(meth_control,meth_test[,i]), cbind(unmeth_control,
                                                                 unmeth_test[,i]), c(s_control,s_test[i]))
  }
  }


  # estimate the abundance of feature
  q0 <- estimateQ((meth_control+unmeth_control), (meth_test+unmeth_test),
                  s_control,s_test)


  # estimate size e
  e_control <- estimateE(meth_control,unmeth_control,s_control,q0)
  if (length(dim(meth_test)[2])==0){
    e_test <- estimateE(meth_test,unmeth_test,s_test,q0)
  }else{
    e_test <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
    for (i in 1:dim(meth_test)[2]){
      e_test[,i] <- estimateE(meth_test[,i],unmeth_test[,i],s_test[i],q0)
    }
  }

  res <- list(p0=p0,p_control=p_control,p_test=p_test,q0=q0,e_control=e_control,
              e_test=e_test)
  return(res)
}


#' The function to calculate all estimate variable for replicate happen on test data.
#'
#' @param meth_control control IP data
#' @param meth_test treated IP data
#' @param unmeth_control control Input data
#' @param unmeth_test treated Input data
#' @param s_control size factor for control data
#' @param s_test size factor for test data
#'
#'
#'
#' @return a list
getBaseEstimate <- function(meth_control,meth_test,unmeth_control,unmeth_test,s_control,s_test){

  # estimate probability of methylation under a condition
  p_control <- estimateP(meth_control, unmeth_control, s_control)
  p_test <- estimateP(meth_test, unmeth_test, s_test)
  p0 <- estimateP(cbind(meth_control,meth_test), cbind(unmeth_control,unmeth_test),
                  c(s_control,s_test))

  # estimate the abundance of feature
  q0 <- estimateQ((meth_control+unmeth_control), (meth_test+unmeth_test),s_control,s_test)

  # estimate size e
  e_control <- estimateE(meth_control,unmeth_control,s_control,q0)
  e_test <- estimateE(meth_test,unmeth_test,s_test,q0)

  res <- list(p0=p0,p_control=p_control,p_test=p_test,q0=q0,e_control=e_control,
              e_test=e_test)
  return(res)
}

#' adjusted p
#'
#' @param x a p vector
#'
#' @return adjust p value
#' @importFrom stats p.adjust
p.adjustfun <-function(x) {p.adjust(x,method="BH")}

#' plot the dispersion fit
#'
#' @param fit_meth fit list of the meth from the locfitW function or
#' the baseFit fucntion
#' @param fit_unmeth fit list of the unmeth from the locfitW function or
#' the baseFit fucntion
#' @param path the path
#'
#' @return a pdf with dispersion fit plot
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics contour par title
plotDispersion <-function(fit_meth,fit_unmeth,path) {
  p <- rep(seq(from=0,to=1,length=101),101)
  q <- rep(seq(from=0,to=100,length=101),101)
  p <- matrix(p,nrow = 101, ncol = 101,byrow=TRUE)
  q <- matrix(q,nrow = 101, ncol = 101,byrow=FALSE)

  p <- p[1:10201]
  q <- q[1:10201]

  w1 <- fittedW(p,q,fit_meth)
  w1 <- matrix(log(w1),nrow = 101, ncol = 101,byrow=FALSE)
  w2 <-  fittedW(p,q,fit_unmeth)
  w2 <- matrix(log(w2),nrow = 101, ncol = 101,byrow=FALSE)


  pdf(path,height=8,width=7)
  #pdf("C:/Users/S41-70/Documents/dispersion.pdf",height=4,width=7)
  par(mfrow=c(2,2))




  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=100,by=0.1),w1)
  title(main="log(w)",sub="meth",ylab="log(q+1)", xlab="p")
  contour(seq(from=0,to=1,by=0.01),seq(from=0,to=100,by=0.1),w2)
  title(main="log(w)",sub="unmeth",ylab="log(q+1)", xlab="p")
  dev.off()
}

#' calculate the p value
#' @param t1 is the total number of detected treated data
#' @param t is the total number of detected data
#' @param n2 is the total number of treated data
#' @param mu2_t mean of t_2
#' @param mu2_c mean of c_2
#' @param size2_t size factor of t_2
#' @param size2_c size factor of c_2
#'
#' @return p-value
#' @importFrom stats dnbinom
quadNBtest <- function(t1,t,n2,mu2_t,mu2_c,size2_t,size2_c){
  nrows <- length(t)
  pval <- rep(1,nrows)


  t2=t-t1



  for (irow in 1:nrows) {

    trip <- n2[irow]

    if (trip<1) {p <- NA} else {

      trip_t2 <- 0:n2[irow]
      trip_c2 <- n2[irow] - trip_t2




      p1 <- dnbinom(x=trip_t2, size=size2_t[irow], mu=mu2_t[irow], log = TRUE)
      if(mu2_t[irow]==0){
        p1 <- rep(0,length(p1))
      }

      p2 <- dnbinom(x=trip_c2, size=size2_c[irow], mu=mu2_c[irow], log = TRUE)
      if(mu2_c[irow]==0){
        p2 <- rep(0,length(p2))
      }
      options(warn=-1)

      # calculate the p value
      p <- p1+p2
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))


      if (is.na(p[(t2[irow]+1)])){
        pv <- NA
      }else{  pv <- sum(p[which(p<=p[(t2[irow]+1)])])}

      pval[irow] <- pv


    }

  }




  res <- data.frame(pval)

  return(res)


}



#' This is the sub function for SIGMRA_similarity_test function
#'
#' @param meth_control detect read counts from false treatment
#' @param meth_test detect read counts from correct treatment
#' @param unmeth_control total read counts minus detection read counts data from false treatment
#' @param unmeth_test total read counts minus detection read counts data from correct treatment
#' @param size.factor TRUE or FALSE
#' @param plot.dispersion whether plot or not
#' @param output.dir the path
#' @param remove.false remove some false positive or not
#' @return a list c("detect proportion treated","detect proportion control","log2 Risk Ratio",
#'  "log2 Odds Ratio","p value","aboundance","adjusted p value")
#' @export
#'
#' @importFrom stats p.adjust
SIGMR_rep_test <-
  function(meth_control,meth_test,unmeth_control,unmeth_test,
           size.factor=NA,
           plot.dispersion=FALSE,
           output.dir = NA,
           remove.false = TRUE) {
    options(warn =-1)
    meth_control <- data.frame(meth_control)
    meth_test <- data.frame(meth_test)
    unmeth_control <- data.frame(unmeth_control)
    unmeth_test <- data.frame(unmeth_test)
    # check
    if( any( ncol(meth_control)!=ncol(unmeth_control), ncol(meth_test)!=ncol(unmeth_test)) ){
      stop( "meth sample and unmeth sample must be the same replicates" )
    }
    if(  ncol(meth_control)<=1 ){
      stop( "number of control sample must be larger than 1" )
    }
    if(  ncol(meth_test)<=1 ){
      stop( "number of test sample must be larger than 1" )
    }
    # estimate
    print("Estimating dispersion for each RNA methylation site, this will take a while ...")




    meth<-cbind(meth_control,meth_test)
    unmeth<-cbind(unmeth_control,unmeth_test)
    if(anyNA(size.factor)){
      s <-sizeFactor(cbind((meth_control+unmeth_control),(meth_test+unmeth_test)))
      s_control <- s[1:length(meth_control[1,])]
      s_test <- s[(length(meth_control[1,])+1):(length(cbind(meth_control,meth_test)[1,]))]
    }else{
      s_control <- size.factor$s_control
      s_test <- size.factor$s_test
    }
    mean <-getBaseEstimate(meth_control,meth_test,unmeth_control,unmeth_test,s_control,s_test)
    p0 <- mean[[1]]
    p1 <- mean[[2]]
    p2 <- mean[[3]]
    q0 <- mean[[4]]
    e1 <- mean[[5]]
    e2 <- mean[[6]]
    res <-baseFit(meth_control,unmeth_control,p1,q0,e1,s_control)
    fit_meth_control<-res[[1]]
    fit_unmeth_control<-res[[2]]
    res <-baseFit(meth_test,unmeth_test,p2,q0,e2,s_test)
    fit_meth_test<-res[[1]]
    fit_unmeth_test<-res[[2]]



    if (is.na(output.dir)) {
      output.dir <- getwd()
    }

    path <- paste(output.dir,"dispersion.pdf",sep = '/')
    if(plot.dispersion){
      plotDispersion(fit_meth_control,fit_unmeth_control,path)
      plotDispersion(fit_meth_test,fit_unmeth_test,path)
    }


    # calculate z
    z_t1 <-calculateZ(q0,p1,s_control,e1)
    z_t2 <-calculateZ(q0,p2,s_test,e2)
    z_c1 <-calculateZ(q0,(1-p1),s_control,e1)
    z_c2 <-calculateZ(q0,(1-p2),s_test,e2)

    # get estimate w
    w_fit_meth_control <-fittedW(p0,q0,fit_meth_control)
    w_fit_meth_test <-fittedW(p0,q0,fit_meth_test)
    w_fit_unmeth_control <-fittedW(p0,q0,fit_unmeth_control)
    w_fit_unmeth_test <-fittedW(p0,q0,fit_unmeth_test)

    # get estimate of upi
    ups_t1 <- pmax(w_fit_meth_control - z_t1, 1e-8)
    ups_t2 <- pmax(w_fit_meth_test - z_t2, 1e-8)
    ups_c1 <- pmax(w_fit_unmeth_control - z_c1, 1e-8)
    ups_c2 <- pmax(w_fit_unmeth_test - z_c2, 1e-8)

    # get all means
    mu_t1 <- (e1*q0*p0)%*%t(as.numeric(s_control))
    mu_t2 <- (e2*q0*p0)%*%t(as.numeric(s_test))
    mu_c1 <- (e1*q0*(1-p0))%*%t(as.numeric(s_control))
    mu_c2 <- (e2*q0*(1-p0))%*%t(as.numeric(s_test))

    # get all variance
    raw_t1 <- (e1%*%t(s_control))^2*ups_t1
    raw_t2 <- (e2%*%t(s_test))^2*ups_t2
    raw_c1 <- (e1%*%t(s_control))^2*ups_c1
    raw_c2 <- (e2%*%t(s_test))^2*ups_c2

    # put mu together
    mu1_t <- rowSums(mu_t1)
    mu2_t <- rowSums(mu_t2)
    mu1_c <- rowSums(mu_c1)
    mu2_c <- rowSums(mu_c2)

    # put size together
    size1_t <- (mu1_t^2)/rowSums(raw_t1)
    size2_t <- (mu2_t^2)/rowSums(raw_t2)
    size1_c <- (mu1_c^2)/rowSums(raw_c1)
    size2_c <- (mu2_c^2)/rowSums(raw_c2)

    # observation together
    t1 <- rowSums(meth_control)
    t2 <- rowSums(meth_test)
    c1 <- rowSums(unmeth_control)
    c2 <- rowSums(unmeth_test)
    t <- t1 + t2
    n1 <- t1 + c1
    n2 <- t2 + c2

    raw <- (rowSums(raw_t1)+rowSums(raw_t2)+rowSums(raw_c1)+rowSums(raw_c2))/4
    # go to test
    res <-quadNBtest(t1,t,n2,mu_t2,mu_c2,
                     size2_t,size2_c)
    if (remove.false==TRUE){
      index <- which(p2<=p1)
      res[index,1] <- 1
    }





    # add fc

    log2.RR <- log2(p2/p1)

    p.treated <- p2
    p.control <- p1

    log2.OR <- log2((rowSums(t(t(meth_test)/s_test))/rowSums(t(t(unmeth_test)/s_test)))/
                      (rowSums(t(t(meth_control)/s_control))/rowSums(t(t(unmeth_control)/s_control))))
    m1 <- rowSums(t(t(meth_control)/s_control))
    m2 <- rowSums(t(t(meth_test)/s_test))
    u1 <- rowSums(t(t(unmeth_control)/s_control))
    u2 <- rowSums(t(t(unmeth_test)/s_test))
    mfc <- log2(m1)-log2(m2)
    ufc <- log2(u1)-log2(u2)

    padj <- p.adjust( res[,1], method="BH" )
    res <- data.frame(p.treated,p.control,log2.RR,log2.OR,res[,1],q0,padj)
    colnames(res) <- c("detect proportion treated","detect proportion control","log2 Risk Ratio",
                       "log2 Odds Ratio","p value","aboundance","adjusted p value")


    return(res)}
