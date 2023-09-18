#' This function serves as the central component for single-cell RNA methylation
#' site analysis by selecting the control group based on the similarity of each
#' test cell. Users can choose the number of control cells, which are selected
#' based on their similarity to the test cells. While it is possible to choose
#' the number of test cells, it is not recommended, and the default is set to 1.
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
#' @param num_test number of cells in the test group (usual be 1 unless you have
#'                          replicates of the test cells)
#' @param num_control number of cells in the control group
#' @param size.factor TRUE or FALSE
#' @param plot.dispersion whether plot or not
#' @param output.dir the path
#' @param remove.false remove some false positive or not
#'
#' @return a list  c("detect proportion treated","detect proportion control","log2 Risk Ratio",
#'  "log2 Odds Ratio","p value","abundance","adjusted p value")
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
#'   method="Pearson", num_test=1,num_control=3,size.factor=NA,plot.dispersion=FALSE,
#'   output.dir = NA,remove.false = TRUE)
#'
#'   # use Euclidean distance
#'   res2 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Euclidean distance",num_test=1,num_control=3,size.factor=NA,
#'   plot.dispersion=FALSE,   output.dir = NA,remove.false = TRUE)
#'
#'  # use expression and Euclidean distance
#'  expression <- matrix (runif(31000,0,1),nrow=1000,ncol=31)
#'  res3<- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'  expression=expression,method="Euclidean distance",num_test=1,num_control=3,
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
#'   method="Pearson",num_test=1,num_control=3,size.factor=NA,plot.dispersion=FALSE,
#'   output.dir = NA,remove.false = TRUE)
#'
#'   # use Euclidean distance
#'   res2 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Euclidean distance",num_test=1,num_control=3,size.factor=NA,
#'   plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
#'
#    # use Pearson with 3 replicates of test data
#'   res3 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Pearson", num_test=3,num_control=3,size.factor=NA,plot.dispersion=
#'   FALSE,output.dir = NA,remove.false = TRUE)
#'
#'   # use Euclidean distance with 3 replicates of test data
#'   res4 <- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'   method="Euclidean distance",num_test=3,num_control=3,size.factor=NA,
#'   plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
#'
#'
#'  # use expression and Euclidean distance
#'  expression <- matrix (runif(37000,0,1),nrow=1000,ncol=37)
#'  res5<- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'  expression=expression,method="Euclidean distance",num_test=1,num_control=3,
#'  size.factor=NA,plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
#'
#'  # use expression and Euclidean distance with 3 replicates of test data
#'  res6<- SIGMR_similarity_test(meth_control,meth_test,unmeth_control,unmeth_test,
#'  expression=expression,method="Euclidean distance",num_test=3,num_control=3,
#'  size.factor=NA,plot.dispersion=FALSE, output.dir = NA,remove.false = TRUE)
SIGMR_similarity_test <-
  function(meth_control,meth_test,unmeth_control,unmeth_test,expression=NA,
           method="Pearson",
           num_test=1,
           num_control,
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
    if(  num_test>ncol(meth_test) ){
      stop( "num_test need smaller or equal to number of cells in the test group" )
    }
    if(  num_control>ncol(meth_control) ){
      stop( "num_control need smaller or equal to number of cells in the control group" )
    }
    if(  num_test<1 ){
      stop( "num_test must be larger or equal to 1" )
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

    if(num_test==1){

      if(method %in% c("Euclidean distance","Manhattan distance","Maximum distance")){

        # calculate the min distance of the test cells in the control group
        similarity_index <- apply(similarity_data,2, function(x){which.minn(x,
                                                        (l_test+l_control))})
        similarity_index <- as.matrix(similarity_index[,(1+l_control):(l_test+l_control)])
        similarity_index_ <- matrix(nrow=l_control,ncol=l_test)
        for (i in 1: l_test){
          similarity_index_[,i] <- similarity_index[which(similarity_index[,
                                            i]%in%c(1:l_control)),i]
        }
        similarity_index <- similarity_index_

        # SIGMRtest

        res <- list()
        for (i in 1:l_test){
          res[[i]] <- SIGMRtest(meth_control[,similarity_index[1:num_control,i]],
          meth_test[,i],
          unmeth_control[,similarity_index[1:num_control,i]],unmeth_test[,i],
          size.factor=NA,plot.dispersion=FALSE,output.dir = NA,
          remove.false = TRUE)
        }
      }else{
        # calculate the min distance of the test cells in the control group
        similarity_index <- apply(similarity_data,2, function(x){which.maxn(x,
                                      (l_test+l_control))})
        similarity_index <- as.matrix(similarity_index[,(1+l_control):(l_test+l_control)])
        similarity_index_ <- matrix(nrow=l_control,ncol=l_test)
        for (i in 1: l_test){
          similarity_index_[,i] <- similarity_index[which(similarity_index[
            ,i]%in%c(1:l_control)),i]
        }
        similarity_index <- similarity_index_

        # SIGMRtest
        res <- list()
        for (i in 1:l_test){
          res[[i]] <- SIGMRtest(meth_control[,similarity_index[1:num_control,i]],
      meth_test[,i],
      unmeth_control[,similarity_index[1:num_control,i]],unmeth_test[,i],
      size.factor=NA,plot.dispersion=FALSE,output.dir = NA,
      remove.false = TRUE)
        }
      }




    }else{

      if(method %in% c("Euclidean distance","Manhattan distance","Maximum distance")){
        # calculate the min distance of the test cells in the control/test group
        similarity_index <- apply(similarity_data,2, function(x){which.minn(x,
                                                        (l_test+l_control))})
        similarity_index <- similarity_index[,(1+l_control):(l_test+l_control)]
        similarity_index_1 <- matrix(nrow=l_control,ncol=l_test)
        similarity_index_2 <- matrix(nrow=l_test,ncol=l_test)
        for (i in 1: l_test){
          similarity_index_1[,i] <- similarity_index[which(similarity_index[,
                                          i]%in%c(1:l_control)),i]
          similarity_index_2[,i] <- similarity_index[which(similarity_index[,i]%in%
                      c((l_control+1):  (l_test+l_control))),i]
        }

        # test
        res <- list()
        for (i in 1:l_test){
          res[[i]] <- SIGMR_rep_test(meth_control[,similarity_index_1[1:num_control,i]],
           meth_test[,(similarity_index_2[1:num_test,i]-dim(meth_control)[2])],
           unmeth_control[,similarity_index_1[1:num_control,i]],
           unmeth_test[,(similarity_index_2[1:num_test,i]-dim(meth_control)[2])],
           size.factor=NA,plot.dispersion=FALSE,output.dir = NA,
           remove.false = TRUE)
        }

      }else{
        # calculate the min distance of the test cells in the control/test group
        similarity_index <- apply(similarity_data,2, function(x){which.maxn(x,
                                                        (l_test+l_control))})
        similarity_index <- similarity_index[,(1+l_control):(l_test+l_control)]
        similarity_index_1 <- matrix(nrow=l_control,ncol=l_test)
        similarity_index_2 <- matrix(nrow=l_test,ncol=l_test)
        for (i in 1: l_test){
          similarity_index_1[,i] <- similarity_index[which(similarity_index[,
                                                    i]%in%c(1:l_control)),i]
          similarity_index_2[,i] <- similarity_index[which(similarity_index[,i]%in%
                           c((l_control+1):  (l_test+l_control))),i]
        }

        # test
        res <- list()
        for (i in 1:l_test){
          res[[i]] <- SIGMR_rep_test(meth_control[,similarity_index_1[1:num_control,i]],
     meth_test[,(similarity_index_2[1:num_test,i]-dim(meth_control)[2])],
     unmeth_control[,similarity_index_1[1:num_control,i]],
     unmeth_test[,(similarity_index_2[1:num_test,i]-dim(meth_control)[2])],
     size.factor=NA,plot.dispersion=FALSE,output.dir = NA,
     remove.false = TRUE)
        }
      }

    }

    res_final <-list()
    for (j in 1:7){
      res_final[[j]] <- data.frame(matrix(nrow=dim(meth_control)[1],ncol=1))
      for (i in 1:l_test){
        res_final[[j]] <- cbind(res_final[[j]],res[[i]][[j]])
      }
      res_final[[j]] <-  res_final[[j]][,-1]
    }

    for (i in 1: length(res_final)){
      res_final[[i]] <- data.frame(res_final[[i]])
      colnames(res_final[[i]]) <- paste0("cell_",(1:dim(res_final[[i]])[2]))
    }
    names(res_final )<- c("detect proportion treated","detect proportion control","log2 Risk Ratio",
                          "log2 Odds Ratio","p value","aboundance","adjusted p value")
    return(res_final)




  }
