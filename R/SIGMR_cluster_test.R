#' This is the main function if you know the cluster of the cells
#' 
#' @description This function serves as the central component for single-cell RNA methylation
#' site analysis by clustering cells based on gene expression or read counts,
#' encompassing both test and control groups. Before calling this function, it
#' is essential to cluster the data. Ideally, each cluster comprises a mix of test
#'  and control cells. Within these clusters, control cell data serves as a
#'  reference, enabling the assessment of test cell behavior. This approach relies
#'  on control cells with similar gene expression patterns as a baseline for
#'  comparison. The method systematically evaluates the behavior of test cells
#'  within their clusters. However, if a cluster lacks control cells, you can
#'  either choose several of the most similar control cells for each test cell
#'  in the cluster (by calling the SIGMR_similarity_test()) or use all control
#'  cells as the background information.
#'
#'
#' @param meth_control detect methylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param meth_test detect methyation read counts from correct treatment
#' (row is gene and column is single cell)
#' @param unmeth_control detect unmethylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param unmeth_test detect unmethylation read counts from false treatment
#' (row is gene and column is single cell)
#' @param expression the expression matrix for similarity, if this variable is miss,
#'                     we use read counts
#' @param cluster cluster index factor vector, so you need cluster the cell before
#'                 put the data into this function
#' @param similarity TRUE or FALSE, if the cluster data do not have control data,
#'                 we use all control data or use similarity to choose some cells.
#' @param method calculate similarity method
#' @param num_control number of cells in the control group (work for similarity)
#' @param size.factor TRUE or FALSE (have size factor or not)
#' @param plot.dispersion whether plot or not
#' @param output.dir the path
#' @param remove.false remove some false positive or not
#' @param withq with the expression abundance or not
#' @param method_dispersion the method to estimate the dispersion ("unbiased_estimate","mle","locfit")
#' @param adjust_1 adjust parameter for RR
#' @param adjust_2 adjust parameter for OR
#' @param one_side p value is one side or both sides
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
#' @export
#'
#'
#' @examples
#'   set.seed(1)
#'   data <- simulateData (test_num=1,control_num=30)
#'
#'   meth_control=data[[1]];meth_test=data[[2]];unmeth_control=data[[3]];unmeth_test=data[[4]]
#'
#'   # simulate the cluster vector
#'   data1 <- cbind(meth_control,meth_test)
#'   num_cluster <- 20
#'   cluster <-factor(sample (1:num_cluster,dim(data1)[2],replace = TRUE))
#'
#'   # cluster with similarity=FALSE
#'   res1 <- SIGMR_cluster_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   cluster=cluster)
#'
#'   # cluster with similarity=TRUE
#'   res2 <- SIGMR_cluster_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   cluster=cluster,similarity=TRUE)
#'
#'   # cluster with expression
#'   expression <- matrix (runif(31000,0,1),nrow=1000,ncol=31)
#'   res3 <- SIGMR_cluster_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   cluster=cluster,expression=expression,similarity=TRUE)
#'
#'
#'
#'
#'   data <- simulateData (test_num=2,control_num=30)
#'
#'   meth_control=data[[1]];meth_test=data[[2]];unmeth_control=data[[3]];unmeth_test=data[[4]]
#'
#'   # simulate the cluster vector
#'   data1 <- cbind(meth_control,meth_test)
#'   num_cluster <- 20
#'   cluster <-factor(sample (1:num_cluster,dim(data1)[2],replace = TRUE))
#'
#'   # cluster with similarity=FALSE
#'   res1 <- SIGMR_cluster_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   cluster=cluster)
#'
#'   # cluster with similarity=TRUE
#'   res2 <- SIGMR_cluster_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   cluster=cluster,similarity=TRUE)
#'
#'   # cluster with expression
#'   expression <- matrix (runif(32000,0,1),nrow=1000,ncol=32)
#'   res3 <- SIGMR_cluster_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   cluster=cluster,
#'   expression=expression,similarity=TRUE)
SIGMR_cluster_test <-
  function(meth_control,meth_test,unmeth_control,unmeth_test,expression=NA,
           cluster,similarity=FALSE,method="Pearson",num_control=3,
           size.factor=NA,
           plot.dispersion=FALSE,
           output.dir = NA,
           remove.false=TRUE,
           withq=TRUE,
           method_dispersion="unbiased_estimate",
           adjust_1=1e-5,
           adjust_2=1e-2,
           one_side=TRUE
  ) {
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
    
    if(  length(expression)!=1){
      if (ncol(expression)!=(ncol(unmeth_control)+ncol(unmeth_test))){
        stop( "col of expression should be the same as the total number of the cells" )
      }
    }
    
    if (length(cluster)!=dim(cbind(meth_control,meth_test))[2]){
      stop( "the length of the cluster should equal to the number of cells" )
      
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
    
    # calculate q
    
    q <- estimateQ((meth_control+unmeth_control), (meth_test+unmeth_test),
                   s_control,s_test)
    
    #SIGMR test
    l <-dim(meth_control)[2]
    
    data1 <- meth_control+unmeth_control
    data2 <- meth_test+unmeth_test
    data <- cbind(data1,data2)
    # test set and control set
    control_set <- c(1:l)
    test_set <- c((l+1):dim(data)[2])
    
    # find the intersection of the cluster and the test set and control set
    
    
    
    
    list_create <- as.list(1:length(levels(cluster)))
    list_function <- function(x){
      index_data <- which(cluster==levels(cluster)[x])
      index_test_data  <- index_data [index_data%in%test_set]
      index_control_data  <- index_data [index_data%in%control_set]
      
      # check whether it contains the control group or not
      if(length(index_control_data )<=1){
        index_null <- x
      }else{
        index_null <-NULL
      }
      return(index_null)
    }
    
    index_null <- unlist(lapply(list_create,list_function))
    
    
    
    
    print(paste0("number of cluster without control group : ",length( index_null )))
    
    # combine the test data without control data
    if(length(index_null)>1){
      levels_new <- levels(cluster)
      levels_new[index_null] <- "null"
      cluster_new <- cluster
      levels(cluster_new)<- levels_new
      
      cluster <- cluster_new
      cluster[which(cluster=="null")[which(cluster=="null")%in%c(
        1:dim(meth_control)[2])]] <- "null_control"
      
    }
    list_create <- as.list(1:length(levels(cluster)))
    list_function1 <- function(x){
      index_data <- which(cluster==levels(cluster)[x])
      index_test_data  <- index_data [index_data%in%test_set]
      index_control_data  <- index_data [index_data%in%control_set]
      
      # check whether it contains the control group or not
      if(length(index_control_data )<=1){
        index_null <- x
      }else{
        index_null <-NULL
      }
      return(index_test_data)
    }
    
    list_function2 <- function(x){
      index_data <- which(cluster==levels(cluster)[x])
      index_test_data  <- index_data [index_data%in%test_set]
      index_control_data  <- index_data [index_data%in%control_set]
      
      # check whether it contains the control group or not
      if(length(index_control_data )<=1){
        index_null <- x
        index_control_data <- control_set
      }else{
        index_null <-NULL
      }
      return(index_control_data)
    }
    
    list_function3 <- function(x){
      index_data <- which(cluster==levels(cluster)[x])
      index_test_data  <- index_data [index_data%in%test_set]
      index_control_data  <- index_data [index_data%in%control_set]
      
      # check whether it contains the control group or not
      if(length(index_control_data )<=1){
        index_null <- x
      }else{
        index_null <-NULL
      }
      return(index_null)
    }
    
    
    index_test_data <- lapply(list_create,list_function1)
    index_control_data <- lapply(list_create,list_function2)
    index_null <- unlist(lapply(list_create,list_function3))
    
    
    
    
    # get the order of the test cells
    list_create <- as.list(1:length(index_test_data))
    list_function <- function(x){
      index_test <- index_test_data [[x]]
    }
    
    index_test <- unlist(lapply(list_create,list_function))-dim(meth_control)[2]
    
    # put into the SIGMRtest function or SIGMR_similarity_test function
    res_list <- list()
    meth <- cbind(meth_control,meth_test)
    unmeth <- cbind(unmeth_control,unmeth_test)
    
    
    list_create <- as.list(1: length(index_test_data))
    
    list_function <- function(x) {
      # if the index_test_data is NA
      if (length(index_test_data[[x]])==0){
        res_list <-NULL
      }else{
        
        # if the test cells do not have control group we use similar cells from the control cells
        if (x %in% index_null & similarity ==TRUE ){
          # if we have the expression we use expression to calculate the similarity
          if (length(expression)==1){
            if(withq ==TRUE){
              q <- q
            } else {q <- FALSE}
            res_list <- SIGMR_similarity_test(meth_control=meth[, index_control_data[[x]]],
                                              meth_test=meth[,index_test_data[[x]]],
                                              unmeth_control=unmeth[,index_control_data[[x]]],unmeth_test=
                                                unmeth[,index_test_data[[x]]],
                                              method=method,num_control=num_control,size.factor,
                                              plot.dispersion, output.dir,
                                              remove.false,q=q, method_dispersion=method_dispersion,adjust_1=adjust_1,
                                              adjust_2=adjust_2,one_side=one_side)
          } else{
            if(withq ==TRUE){
              q <- q
            } else {q <- FALSE}
            # if we do not have the expression we use read counts to calculate the similarity
            res_list<- SIGMR_similarity_test(meth_control=meth[,index_control_data[[x]]]
                                             ,meth_test=meth[,index_test_data[[x]]],
                                             unmeth_control=unmeth[,index_control_data[[x]]],unmeth_test=
                                               unmeth[,index_test_data[[x]]],
                                             expression=expression[,c(index_control_data[[x]],index_test_data[[x]])],
                                             method=method,num_control=num_control,size.factor,
                                             plot.dispersion, output.dir,
                                             remove.false,q=q,method_dispersion=method_dispersion,adjust_1=adjust_1,
                                             adjust_2=adjust_2 ,one_side=one_side)
          }
          
        }else{
          if(withq ==TRUE){
            q <- q
          } else {q <- FALSE}
          res_list <- SIGMRtest(meth_control=meth[,index_control_data[[x]]]
                                ,meth_test=meth[,index_test_data[[x]]],
                                unmeth_control=unmeth[,index_control_data[[x]]],
                                unmeth_test=unmeth[,index_test_data[[x]]],
                                size.factor,     plot.dispersion, output.dir,
                                remove.false,q=q ,method_dispersion=method_dispersion,adjust_1=adjust_1,
                                adjust_2=adjust_2,one_side=one_side)
        }
        
      }
      return(res_list)
    }
    
    res_list <-  lapply(list_create,list_function)
    
    
    
    
    list_create <- as.list(c(1,3,4,5,7,8))
    
    list_function <- function(x) {
      
      list_create_ <- as.list(1: length(index_test_data))
      
      list_function_ <- function(y){
        res <- res_list[[y]][[x]]
      }
      res <- lapply(list_create_,list_function_)
      res <- as.data.frame(res[which(lapply(res,is.null)==FALSE)])
      res <- res[, order(index_test)]
      
    }
    
    res1 <- lapply(list_create,list_function)
    
    
    
    
    
    
    list_create <- as.list(c(2,6))
    
    list_function <- function(x) {
      
      list_create_ <- as.list(1: length(index_test_data))
      
      list_function_ <- function(y){
        
        if (is.null(res_list[[y]])){
          res <- NULL
        }else{
          res <- matrix(rep(res_list[[y]][[x]],(dim(res_list[[y]][[1]])[2])),nrow=
                          (dim(res_list[[y]][[1]])[1]),
                        ncol=(dim(res_list[[y]][[1]])[2]))
        }
        
      }
      res <- lapply(list_create_,list_function_)
      res <- as.data.frame(res[which(lapply(res,is.null)==FALSE)])
      res <- res[, order(index_test)]
      
    }
    
    res2 <- lapply(list_create,list_function)
    
    
    res <- c(res1,res2)[order(c(1,3,4,5,7,8,2,6))]
    
    
    
    list_create <- as.list(1: length(res))
    list_function <- function(x) {
      if(is.null(dim(res[[x]])[2])){
      }else{
        colnames(res[[x]]) <- paste0("cell_",(1:dim(res[[x]])[2]))
      }
      return(res[[x]])
    }
    
    
    res_final <- lapply(list_create,list_function)
    
    
    names(res_final )<- c("detect proportion treated","detect proportion control",
                          "log2 Risk Ratio","log2 Odds Ratio","p value","aboundance",
                          "adjusted p value","TCR")
    
    return(res_final )
    
    
    
  }
