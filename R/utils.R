## This is the document that contains the sub function for the main function

#' @importFrom stats median 
#' @import locfit
#' @importFrom stats predict optimize optimise
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics contour par title
#' @importFrom stats dnbinom
## get size factor function.

rowMeans_function <- function(data){
  l <- apply(data,1,function(x)length(which(x!=0)))
  mean_value <- rowSums(data)/l
  mean_value [which(l==0)] <- 0
  mean_value
}
sizeFactor <- function(data) {
  # make the elements no smaller than 0
  data <- as.matrix(data)
  
  # first log
  log_data <- log(data)
  log_data[is.infinite(log_data)] <- NA
  log_mean <- rowMeans(log_data)
  log_s <- log_data-log_mean
  
  # then exp
  s_size <- exp(apply(log_s,2,function(x)median(x,na.rm=TRUE)))
  return(s_size)
}



## estimate P function.
estimateP <- function(meth, unmeth, size) {
  p <- rowSums(t(t(meth)/size))/rowSums(t(t(unmeth+meth)/size))
  return(p)
}

## estimate Q function.

estimateQ <- function(total_control,total_test,s_control,s_test){
  temp <- t(t(cbind(total_control,total_test))/c(s_control,s_test))
  q <- rowMeans_function(temp)
  return(q)
}

## estimate E function.

estimateE <- function(meth,unmeth,size,q){
  
  e <- rowMeans_function(t(t(meth+unmeth)/size))/q
  return(e)
}

## The function to get estimated v

calculate_v <- function(data,size,e){
  
  q <-  t(t(data)/size)/e
  q <- q-rowMeans_function(q)
  v <- rowSums(q^2)/(length(size)-1)
  v <- as.matrix(v)
  return(v)
}


## The function to get estimated Z

calculateZ <- function(q,p,size,e){
  temp <- p*q/length(size)
  
  z <- rowSums(1/e%*%t(size))*temp
  
  
  return(z)
}


## The function to local fit the v

locfit_v <- function(p,q,v) {
  l <- log(q+1)
  intercept<-c(rep(1,length(q)))
  data<-data.frame(cbind(p,l,intercept,v))
  fit=locfit(v~lp(p,l,intercept),data=data,family="gamma")
  return(fit)
}




## The function to get final estimated Z

fitted_v <- function(p,q,fit){
  l <- log(q+1)
  intercept<-rep(1,length(q))
  data=data.frame(cbind(p,l,intercept))
  index <- which(is.na(p))
  if(length(index)==0){
    v_fit <- predict(fit,data)
  }else{
    v_fit <- rep(NA,length(q))
    v_fit[-index] <- predict(fit,data[-index,])
  }
  
  return(v_fit)
}



## The function to get the fit of W

baseFit <- function(meth,unmeth,p,q,e,size){
  options(warn =-1)
  # get estimated w
  v_meth <-calculate_v(meth,size,e)
  v_unmeth <-calculate_v(unmeth,size,e)
  
  # locfit
  fit_meth <- locfit_v(p,q,v_meth)
  fit_unmeth <- locfit_v((1-p),q,v_unmeth)
  
  res <- list(fit_meth,fit_unmeth)
  return(res)
}

## The function to calculate all estimate variable.

getEstimate <- function(meth_control,meth_test,unmeth_control,unmeth_test,s_control,s_test){
  
  # estimate probability of methylation under a condition
  p_control <- estimateP(meth_control, unmeth_control, s_control)
  if (length(dim(meth_test)[2])==0){
    
    p_test <- estimateP(meth_test, unmeth_test, s_test)
    p0 <- estimateP(cbind(meth_control,meth_test), cbind(unmeth_control,
                                                         unmeth_test), c(s_control,s_test))
    
  }else{p_test <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
  p0 <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],length(s_test)),
                                             as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
  
  list_function <- function(x) {
    estimateP( x[2:(x[1]+1)], x[(2+x[1]):(2*x[1]+1)], x[2*x[1]+2])
  }
  
  p_test <- as.data.frame(lapply(list_create,list_function))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],length(s_test)),
                                             as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
  
  list_function <- function(x) {
    estimateP( cbind(meth_control,x[2:(x[1]+1)]), cbind(unmeth_control,x[(2+x[1]):
                                                                           (2*x[1]+1)]), c(s_control,x[2*x[1]+2]))
  }
  
  p0 <- as.data.frame(lapply(list_create,list_function))
  
  
  }
  
  
  # estimate the abundance of feature
  q0 <- estimateQ((meth_control+unmeth_control), (meth_test+unmeth_test),
                  s_control,s_test)
  
  
  # estimate size e
  e_control <- estimateE(meth_control,unmeth_control,s_control,q0)
  if (length(dim(meth_test)[2])==0){
    e_test <- estimateE(meth_test,unmeth_test,s_test,q0)
  }else{
    
    list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],length(s_test)),
                                               as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
    
    list_function <- function(x) {
      estimateE( x[2:(x[1]+1)], x[(2+x[1]):(2*x[1]+1)], x[2*x[1]+2],q0)
    }
    
    e_test <- as.data.frame(lapply(list_create,list_function))
  }
  
  res <- list(p0=p0,p_control=p_control,p_test=p_test,q0=q0,e_control=e_control,
              e_test=e_test)
  return(res)
}



getEstimate_withq <- function(meth_control,meth_test,unmeth_control,unmeth_test,
                              s_control,s_test,q){
  
  # estimate probability of methylation under a condition
  p_control <- estimateP(meth_control, unmeth_control, s_control)
  if (length(dim(meth_test)[2])==0){
    
    p_test <- estimateP(meth_test, unmeth_test, s_test)
    p0 <- estimateP(cbind(meth_control,meth_test), cbind(unmeth_control,
                                                         unmeth_test), c(s_control,s_test))
    
  }else{p_test <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
  p0 <- data.frame(matrix(nrow=dim(meth_test)[1],ncol=dim(meth_test)[2]))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],length(s_test)),
                                             as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
  
  list_function <- function(x) {
    estimateP( x[2:(x[1]+1)], x[(2+x[1]):(2*x[1]+1)], x[2*x[1]+2])
  }
  
  p_test <- as.data.frame(lapply(list_create,list_function))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],length(s_test)),
                                             as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
  
  list_function <- function(x) {
    estimateP( cbind(meth_control,x[2:(x[1]+1)]), cbind(unmeth_control,x[(2+x[1]):
                                                                           (2*x[1]+1)]), c(s_control,x[2*x[1]+2]))
  }
  
  p0 <- as.data.frame(lapply(list_create,list_function))
  
  
  }
  
  
  # estimate size e
  e_control <- estimateE(meth_control,unmeth_control,s_control,q)
  if (length(dim(meth_test)[2])==0){
    e_test <- estimateE(meth_test,unmeth_test,s_test,q)
  }else{
    
    list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],length(s_test)),
                                               as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
    
    list_function <- function(x) {
      estimateE( x[2:(x[1]+1)], x[(2+x[1]):(2*x[1]+1)], x[2*x[1]+2],q)
    }
    
    e_test <- as.data.frame(lapply(list_create,list_function))
  }
  
  res <- list(p0=p0,p_control=p_control,p_test=p_test,e_control=e_control,
              e_test=e_test)
  return(res)
}




## adjusted p
p.adjustfun <-function(x) {p.adjust(x,method="BH")}

## plot the dispersion fit

plotDispersion <-function(fit_meth,fit_unmeth,path) {
  
  
  pdf(path,height=8,width=7)
  #pdf("C:/Users/S41-70/Documents/dispersion.pdf",height=4,width=7)
  if(inherits(fit_meth,"locfit")){
    plot(fit_meth)
  }else{
    print("locfit for fit_meth can not plot")
  }
  if(inherits(fit_unmeth,"locfit")){
    plot(fit_unmeth)
  }else{
    print("locfit for fit_unmeth can not plot")
  }
  dev.off()
}

## calculate the p value

quadNBtest <- function(t2,n2,mu2_t,mu2_c,size2_t,size2_c,one_side=TRUE){
  
  
  
  
  list_create <- as.list(as.data.frame(t(cbind(t2,n2,mu2_t,mu2_c,size2_t,size2_c))))
  
  
  list_function <- function(x){
    
    trip <- x[2]
    if ( trip < 1) {
      pv <- NA
    }else {
      
      trip_t2 <- 0:x[2]
      trip_c2 <- x[2] -  trip_t2
      
      
      
      
      p1 <- dnbinom(x=trip_t2, size=x[5], mu=x[3], log = TRUE)
      if(x[3]==0){
        p1 <- rep(0,length(p1))
      }
      
      p2 <- dnbinom(x=trip_c2, size=x[6], mu=x[4], log = TRUE)
      if(x[4]==0){
        p2 <- rep(0,length(p2))
      }
      options(warn=-1)
      
      # calculate the p value
      p <- p1+p2
      p <- p-max(p)
      p <- exp(p)/sum(exp(p))
      
      
      if (is.na(p[x[1]+1])){
        pv <- NA
      }else{  
        if(one_side==TRUE){
          pv <- sum(p[which(trip_t2>=x[1]+1)])
        }else{
          pv <- sum(p[which(p<=p[(x[1]+1)])])
        }
        
        
        
      }
      
      
    }
    return(pv)
    
  }
  
  
  
  res <- unlist(lapply(list_create,list_function)) 
  
  return(res)
  
  
}






Estimate_dispersion <- function(data_meth,data_unmeth,mu_meth,mu_unmeth){
  
  if( any( ncol(data_meth)!=ncol(mu_meth), ncol(data_unmeth)!=ncol(mu_unmeth)) ){
    stop( "data parameters and mu sparameters need to have same number of column" )
  }
  if( any( nrow(data_meth)!=nrow(data_unmeth), nrow(data_meth)!=nrow(mu_meth), 
           nrow(mu_meth)!=nrow(mu_unmeth)) ){
    stop( "four samples must have same number of sites" )
  }
  
  dim1 <- dim(data_meth)[1]
  dim2 <- dim(data_meth)[2]
  dim3 <- dim(data_unmeth)[2]
  
  
  meth_list <- as.list(as.data.frame(rbind(t(data_meth),t(mu_meth))))
  unmeth_list <- as.list(as.data.frame(rbind(t(data_unmeth),t(mu_unmeth))))
  
  
  
  
  function_estimate <- function(df)  {
    function_loglikelihood <- function(alpha){
      sum(dnbinom(df[1:dim2],mu=df[(1+dim2):(2*dim2)],size=alpha,log = TRUE),na.rm=TRUE)}
    
    op <- optimise(function_loglikelihood,c(1e-10,1e5),maximum = TRUE)
    op$maximum
  }
  
  res <- list()
  res[[1]] <- unlist(lapply(meth_list, function_estimate))
  
  res[[2]] <- unlist(lapply(unmeth_list, function_estimate))
  
  
  return(res)
  
}
