#' The function to get the three quantification levels
#'
#' @param meth_control detect methylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param meth_test detect methyation read counts from correct treatment
#' (row is gene and column is single cell)
#' @param unmeth_control detect unmethylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param unmeth_test detect unmethylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param size.factor NA or a data.frame with one column is control size factor and test size factor
#' @param q expression aboundance
#' @param adjust_1 adjust parameter for RR
#' @param adjust_2 adjust parameter for OR
#'
#' @return a list c("methylation level treated","methylation level control","log2 Risk Ratio",
#' "log2 Odds Ratio","TCR")
#'
#'  1.The first part of the list is the meth proportion data frame (row is gene
#'  and column is single cell) in the test cells.
#'  2.The second part of the list is the mean meth proportion vector in the control cells.
#'  3.The third part of the list is the log2 risk ratio data frame for the test cells.
#'  4.The forth part of the list is the log2 odds ratio data frame for the test cells.
#'  5.The fifth part of the list is the TCR for the test cells.
#'  
#'  For NA in the return list:
#'  
#'  1. NA in the treated methylation level, it is caused by that the test cells have no reads 
#'  in both methylation reads and un methylation reads
#'  2. NA in the control methylation level, it is caused by that the all control cells have 
#'  no reads in both methylation reads and un methylation reads
#'  3. NA in the risk ratio, the methylation proportion is NA on the test cell or control cells group.
#'  4. NA in the TCR, is similar to the risk ratio.
#' @export
#'
#' @examples 
#'  data <- simulateData (test_num=2,control_num=30)
#'  res <- quantification_level(data[[1]],data[[2]],data[[3]],data[[4]])
quantification_level <- function(meth_control,meth_test,unmeth_control,unmeth_test,
                                 size.factor=NA,q=NA,adjust_1=1e-5,adjust_2=1e-2){
  options(warn =-1)
  meth_control <- data.frame(meth_control)
  meth_test <- data.frame(meth_test)
  unmeth_control <- data.frame(unmeth_control)
  unmeth_test <- data.frame(unmeth_test)
  # check
  if( any( ncol(meth_control)!=ncol(unmeth_control), ncol(meth_test)!=ncol(unmeth_test)) ){
    stop( "meth sample and unmeth sample must be the same replicates" )
  }
  if( any( nrow(meth_control)!=nrow(meth_test), nrow(meth_test)!=nrow(unmeth_test), 
           nrow(meth_test)!=nrow(unmeth_control)) ){
    stop( "mall four samples must have same number of sites" )
  }
  
  if(  ncol(meth_control)<=1 ){
    stop( "number of control sample must be larger than 1" )
  }
  
  
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
  
  if (length(q)==1){
    mean <-getEstimate(meth_control,meth_test,unmeth_control,unmeth_test,s_control,s_test)
    p0 <- mean[[1]]
    p_control <- mean[[2]]
    p_test <- mean[[3]]
    q0 <- mean[[4]]
    e_control <- mean[[5]]
    e_test <- mean[[6]]
  }else{
    mean <-getEstimate_withq(meth_control,meth_test,unmeth_control,unmeth_test,
                             s_control,s_test,q)
    p0 <- mean[[1]]
    p_control <- mean[[2]]
    p_test <- mean[[3]]
    q0 <- q
    e_control <- mean[[4]]
    e_test <- mean[[5]]
  }
  
  
  
  p.treated <- p_test
  p.control <- p_control
  
  # add RR 
  log2.RR <- log2((p_test+adjust_1)/(p_control+adjust_1))
  
  
  list_create <- as.list(as.data.frame(rbind(rep(dim(meth_test)[1],dim(meth_test)[2]),
                                             as.matrix(meth_test),as.matrix(unmeth_test),s_test)))
  
  
  list_function <- function(x){
    numerator <- ((x[2:(1+x[1])]/x[2*x[1]+2])*
                    (rowMeans(t(t(unmeth_control)/s_control))))/(q0)^2
    denominator <- ((rowMeans(t(t(meth_control)/s_control)))* 
                      (x[(2+x[1]):(1+2*x[1])]/x[2*x[1]+2]))/(q0)^2
    adjust <- max(numerator[which(denominator==0)],na.rm=TRUE)/2^10
    
    log2.OR <- log2((numerator+adjust*adjust_2)/(denominator+adjust*adjust_2))
    log2.OR
  }
  
  
  log2.OR <- as.data.frame(lapply(list_create, list_function))
  
  list_create <- as.list(1:dim(log2.OR)[2])
  
  list_function <- function(x){
    log2.OR[which(is.na(log2.RR[,x])),x] <- NA
    log2.OR[,x]
  }
  
  log2.OR <- as.data.frame(lapply(list_create, list_function))
  
  
  TCR <- p.treated-p.control
  
  TCR <- data.frame(TCR)
  
  p.treated <- data.frame(p.treated)
  colnames(p.treated) <- paste0("cell_",(1:dim(p.treated)[2]))
  
  log2.RR <- data.frame(log2.RR)
  colnames(log2.RR) <- paste0("cell_",(1:dim(log2.RR)[2]))
  
  log2.OR <- data.frame(log2.OR)
  colnames(log2.OR) <- paste0("cell_",(1:dim(log2.OR)[2]))
  
  TCR <- data.frame(TCR)
  colnames(TCR) <- paste0("cell_",(1:dim(TCR)[2]))
  
  res <-list(p.treated,p.control,log2.RR,log2.OR,TCR)
  
  
  
  names(res) <-c("methylation level treated","methylation level control","log2 Risk Ratio",
                 "log2 Odds Ratio","TCR")
  return(res)
  
  
}