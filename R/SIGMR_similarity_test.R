#' This is the function work for the SIGMR_cluster_test
#' 
#' @description This function serves as the central component for single-cell RNA methylation
#' site analysis by selecting the control group based on the similarity of each
#' test cell. Users can choose the number of control cells, which are selected
#' based on their similarity to the test cells. 
#'
#' @param meth_control detect methylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param meth_test detect methyation read counts from correct treatment
#' (row is gene and column is single cell)
#' @param unmeth_control detect unmethylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param unmeth_test detect unmethylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param expression the expression matrix for compare the similarity, if this
#'                              variable is miss, we use read counts
#' @param method calculate similarity method c("Euclidean distance","Manhattan distance",
#'         "Maximum distance","Pearson","Spearman")
#' @param num_control number of cells in the control group
#' @param size.factor TRUE or FALSE
#' @param plot.dispersion whether plot or not
#' @param output.dir the path
#' @param remove.false remove some false positive or not
#' @param q expression aboundance 
#' @param method_dispersion the method to estimate the dispersion ("unbiased_estimate","mle","locfit")
#' @param adjust_1 adjust parameter for RR
#' @param adjust_2 adjust parameter for OR
#' @param one_side p value is one side or both sides
#' 
#' 
#' @return a list c("methylation level treated","methylation level control","log2 Risk Ratio",
#' "log2 Odds Ratio","p value","aboundance","adjusted p value","TCR")
#'
#'  1.The first part of the list is the meth proportion data frame (row is gene
#'  and column is single cell) in the test cells.
#'  2.The second part of the list is the mean meth proportion data frame (row is gene
#'  and column is single cell) in the control cells.
#'  3.The third part of the list is the log2 risk ratio data frame for the test cells.
#'  4.The forth part of the list is the log2 odds ratio data frame for the test cells.
#'  5.The fifth part of the list is the p value data frame for the test cells.
#'  6.The sixth part of the list is the estimated gene abundance data frame.
#'  7.The seventh part of the list is the adjusted p value data frame for the test cells.
#'  8.The eighth part of the list is the TCR for the test cells.
#'  
#'  For NA in the return list:
#'  
#'  1. NA in the treated methylation level, it is caused by that the test cells have no reads 
#'  in both methylation reads and un methylation reads
#'  2. NA in the control methylation level, it is caused by that the all control cells have 
#'  no reads in both methylation reads and un methylation reads
#'  3. NA in the risk ratio, the methylation proportion is NA on the test cell or control cells group.
#'  4. NA in the TCR, is similar to the risk ratio.
#'  5. NA in the p value, it is only caused by the that control methylation level is NA and the 
#'  treated methylation level is not NA.
#'  
#' @export
#'
#' @importFrom doBy which.minn which.maxn
#' @importFrom stats cor dist sd
#'
#' @examples
#'   data <- simulateData (test_num=1,control_num=30)
#'
#'   meth_control=data[[1]];meth_test=data[[2]];unmeth_control=data[[3]];unmeth_test=data[[4]]
#'
#'   # use Pearson
#'   res1 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Pearson",num_control=3,size.factor=NA,plot.dispersion=FALSE,
#'   output.dir = NA,remove.false = TRUE)
#'
#'   # use Euclidean distance
#'   res2 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Euclidean distance",num_control=3,size.factor=NA,
#'   plot.dispersion=FALSE,   output.dir = NA,remove.false = TRUE)
#'
#'  # use expression and Euclidean distance
#'  expression <- matrix (runif(31000,0,1),nrow=1000,ncol=31)
#'  res3<- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'  expression=expression,method="Euclidean distance",num_control=3,
#'  size.factor=NA,plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
#'
#'
#'
#'   data <- simulateData (test_num=7,control_num=30)
#'
#'   meth_control=data[[1]];meth_test=data[[2]];unmeth_control=data[[3]];unmeth_test=data[[4]]
#'
#'   # use Pearson
#'   res1 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Pearson",num_control=3,size.factor=NA,plot.dispersion=FALSE,
#'   output.dir = NA,remove.false = TRUE)
#'
#'   # use Euclidean distance
#'   res2 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Euclidean distance",num_control=3,size.factor=NA,
#'   plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
#'
#'
#'
#'  # use expression and Euclidean distance
#'  expression <- matrix (runif(37000,0,1),nrow=1000,ncol=37)
#'  res5<- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'  expression=expression,method="Euclidean distance",num_control=3,
#'  size.factor=NA,plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
#'
SIGMR_similarity_test <-
  function(meth_control,meth_test,unmeth_control,unmeth_test,expression=NA,
           method="Pearson",
           num_control,
           size.factor=NA,
           plot.dispersion=FALSE,
           output.dir = NA,
           remove.false = TRUE,
           q=NA,
           method_dispersion="unbiased_estimate",
           adjust_1=1e-5,
           adjust_2=1e-2,
           one_side=TRUE) {
    options(warn =-1)
    meth_control <- data.frame(meth_control)
    meth_test <- data.frame(meth_test)
    unmeth_control <- data.frame(unmeth_control)
    unmeth_test <- data.frame(unmeth_test)
    # check
    if( any( ncol(meth_control)!=ncol(unmeth_control), ncol(meth_test)!=ncol(
      unmeth_test)) ){
      stop( "meth sample and unmeth sample must be the same replicates" )
    }
    if( any( nrow(meth_control)!=nrow(meth_test), nrow(meth_test)!=nrow(unmeth_test), 
             nrow(meth_test)!=nrow(unmeth_control)) ){
      stop( "mall four samples must have same number of sites" )
    }
    if(  ncol(meth_control)<=1 ){
      stop( "number of control sample must be larger than 1" )
    }
    if(  num_control<=1 ){
      stop( "num_control must be larger than 1" )
    }
    
    if(  num_control>ncol(meth_control) ){
      stop( "num_control need smaller or equal to number of cells in the control group" )
    }
    
    
    if(  length(expression)!=1){
      if (ncol(expression)!=(ncol(unmeth_control)+ncol(unmeth_test))){
        stop( "col of expression should be the same as the total number of the cells" )
      }
    }
    
    
    
    
    
    # calculate log normalized data and scale by z score
    
    data1 <- meth_control+unmeth_control
    data2 <- meth_test+unmeth_test
    l_test <- dim(data2)[2]
    l_control <- dim(data1)[2]
    if (length(expression)==1){
      data <- cbind(data1,data2)
      
      data <- apply(data,2,function(x) log2(x/sum(x)*10000+1))
      data <- apply(data,2,function(x) (x-mean(x))/sd(x))
    }else{
      data <- expression
    }
    
    
    
    #data1 <- data[,1:l_control ]
    #data2 <- data[,(l_control+1): (l_control+l_test)]
    
    # calculate similarity
    similarity_data <- matrix(nrow=(l_test+l_control),ncol=(l_test+l_control))
    if (method=="Euclidean distance"){
      
      similarity_data <-  as.matrix(dist(t(as.matrix(data)), method = "euclidean",
                                         diag = TRUE, upper = TRUE))
      
    } else if (method=="Manhattan distance"){
      
      similarity_data <- as.matrix(dist(t(as.matrix(data)), method = "manhattan",
                                        diag = TRUE, upper = TRUE))
      
    } else if (method=="Maximum distance"){
      
      similarity_data <- as.matrix(dist(t(as.matrix(data)), method = "maximum",
                                        diag = TRUE, upper = TRUE))
      
    } else if (method=="Pearson"){
      
      similarity_data <- cor(data, method = "pearson")
      
    } else if (method=="Spearman"){
      
      
      similarity_data <- cor(data, method = "spearman")
      
    }
    
    # test
    
    
    
    if(method %in% c("Euclidean distance","Manhattan distance","Maximum distance")){
      
      # calculate the min distance of the test cells in the control group
      similarity_index <- apply(similarity_data,2, function(x){which.minn(x,
                                                                          (l_test+l_control))})
      similarity_index <- as.matrix(similarity_index[,(1+l_control):(l_test+l_control)])
      
      list_create <- as.list(1: l_test)
      
      list_function <- function(x){
        similarity_index_ <- similarity_index[which(similarity_index[,
                                                                     x]%in%c(1:l_control)),x]
        return(similarity_index_)
      }
      
      similarity_index <- as.data.frame(lapply(list_create,list_function))
      
      # SIGMRtest
      
      list_create <- as.list(1: l_test)
      
      list_function <- function(x){
        res <- SIGMRtest(meth_control[,similarity_index[1:num_control,x]],
                         meth_test[,x],
                         unmeth_control[,similarity_index[1:num_control,x]],unmeth_test[,x],
                         size.factor=NA,plot.dispersion=FALSE,output.dir = NA,
                         remove.false = TRUE,q=q,
                         method_dispersion=method_dispersion,adjust_1=adjust_1,
                         adjust_2=adjust_2,one_side=one_side)
        return(res)
      }
      
      res <- lapply(list_create,list_function)
      
    }else{
      # calculate the max similarity of the test cells in the control group
      similarity_index <- apply(similarity_data,2, function(x){which.maxn(x,
                                                                          (l_test+l_control))})
      similarity_index <- as.matrix(similarity_index[,(1+l_control):(l_test+l_control)])
      list_create <- as.list(1: l_test)
      
      list_function <- function(x){
        similarity_index_ <- similarity_index[which(similarity_index[,
                                                                     x]%in%c(1:l_control)),x]
        return(similarity_index_)
      }
      
      similarity_index <- as.data.frame(lapply(list_create,list_function))
      
      # SIGMRtest
      list_create <- as.list(1: l_test)
      
      list_function <- function(x){
        res <- SIGMRtest(meth_control[,similarity_index[1:num_control,x]],
                         meth_test[,x],
                         unmeth_control[,similarity_index[1:num_control,x]],unmeth_test[,x],
                         size.factor=NA,plot.dispersion=FALSE,output.dir = NA,
                         remove.false = TRUE,q=q,
                         method_dispersion=method_dispersion,adjust_1=adjust_1,
                         adjust_2=adjust_2,one_side=one_side)
        return(res)
      }
      
      res <- lapply(list_create,list_function)
      
    }
    
    
    
    
    list_create <- as.list(1: 8)
    list_function <- function(x){
      list_create_ <- as.list(1: l_test)
      list_function_ <- function(y){
        res[[y]][[x]]
      }
      
      res <- as.data.frame(lapply(list_create_,list_function_))
    }
    
    res_final <- lapply(list_create,list_function)
    
    list_create <- as.list(1: length(res_final))
    list_function <- function(x){
      res_final[[x]] <- data.frame(res_final[[x]])
      
      
      if(is.null(dim(res_final[[x]])[2])){
      }else{
        colnames(res_final[[x]]) <- paste0("cell_",(1:dim(res_final[[x]])[2]))
      }
      return(res_final[[x]])
    }
    
    res_final <- lapply(list_create,list_function)
    
    
    names(res_final )<- c("detect proportion treated","detect proportion control",
                          "log2 Risk Ratio","log2 Odds Ratio","p value","aboundance",
                          "adjusted p value","TCR")
    return(res_final)
    
    
    
    
  }
