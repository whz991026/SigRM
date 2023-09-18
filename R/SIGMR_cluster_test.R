#' This function serves as the central component for single-cell RNA methylation
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
#' @param num_test number of cells in the test group
#' @param num_control number of cells in the control group (work for similarity)
#' @param size.factor TRUE or FALSE (have size factor or not)
#' @param plot.dispersion whether plot or not
#' @param output.dir the path
#' @param remove.false remove some false positive or not
#'
#' @return a list  c("detect proportion treated","detect proportion control","log2 Risk Ratio",
#' "log2 Odds Ratio","p value","abundance","adjusted p value")
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
           cluster,similarity=FALSE,method="Pearson",
           num_test=1,num_control=3,
           size.factor=NA,
           plot.dispersion=FALSE,
           output.dir = NA,
           remove.false=TRUE
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

    # cluster



    l <-dim(meth_control)[2]

    data1 <- meth_control+unmeth_control
    data2 <- meth_test+unmeth_test
    data <- cbind(data1,data2)


    #SIGMR test

    # test set and control set
    control_set <- c(1:l)
    test_set <- c((l+1):dim(data)[2])

    # find the intersection of the cluster and the test set and control set
    index_data <- list()
    index_test_data <- list()
    index_control_data <- list()
    index_null <- c()
    for (i in 1:length(levels(cluster))){

      index_data[[i]] <- which(cluster==levels(cluster)[i])
      index_test_data [[i]] <- index_data [[i]][index_data[[i]]%in%test_set]
      index_control_data [[i]] <- index_data [[i]][index_data[[i]]%in%control_set]

      # check whether it contains the control group or not
      if(length(index_control_data [[i]])<=1){
        index_null <- c(index_null,i)
      }

    }
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
    index_data <- list()
    index_test_data <- list()
    index_control_data <- list()
    index_null <- c()
    for (i in 1:length(levels(cluster))){

      index_data[[i]] <- which(cluster==levels(cluster)[i])
      index_test_data [[i]] <- index_data [[i]][index_data[[i]]%in%test_set]
      index_control_data [[i]] <- index_data [[i]][index_data[[i]]%in%control_set]

      # check whether it contains the control group or not
      if(length(index_control_data [[i]])<=1){



        print("some test data have not enough control group,use the whole group")
        index_null <- c(index_null,i)
        print(paste0("the index is ",i))
        index_control_data [[i]] <- control_set



      }


    }



    # get the order of the test cells
    index_test <-c()
    for (i in 1: length(index_test_data)){
      index_test_ <- index_test_data [[i]]
      index_test <-c(index_test,index_test_)
    }
    index_test <- index_test -dim(meth_control)[2]

    # put into the SIGMRtest function or SIGMR_similarity_test function
    res_list <- list()
    meth <- cbind(meth_control,meth_test)
    unmeth <- cbind(unmeth_control,unmeth_test)



    for ( i in 1: length(index_test_data)){
      # if the index_test_data is NA
      if (length(index_test_data[[i]])==0){
        res_list[[i]] <-list()
      }else{

        # if the test cells do not have control group we use similar cells from the control cells
        if (i %in% index_null & similarity ==TRUE ){
          # if we have the expression we use expression to calculate the similarity
          if (length(expression)==1){
            res_list[[i]] <- SIGMR_similarity_test(meth_control=meth[,
                           index_control_data[[i]]],
         meth_test=meth[,index_test_data[[i]]],
         unmeth_control=unmeth[,index_control_data[[i]]],unmeth_test=
           unmeth[,index_test_data[[i]]],
         method=method, num_test=num_test,num_control=num_control,size.factor,
         plot.dispersion, output.dir,
         remove.false )
          } else{
            # if we do not have the expression we use read counts to calculate the similarity
            res_list[[i]] <- SIGMR_similarity_test(meth_control=meth[,index_control_data[[i]]]
           ,meth_test=meth[,index_test_data[[i]]],
           unmeth_control=unmeth[,index_control_data[[i]]],unmeth_test=
             unmeth[,index_test_data[[i]]],
           expression=expression[,c(index_control_data[[i]],index_test_data[[i]])],
           method=method, num_test=num_test,num_control=num_control,size.factor,
           plot.dispersion, output.dir,
           remove.false )
          }

        }else{
          res_list[[i]] <- SIGMRtest(meth_control=meth[,index_control_data[[i]]]
                                     ,meth_test=meth[,index_test_data[[i]]],
                                     unmeth_control=unmeth[,index_control_data[[i]]],
                                     unmeth_test=unmeth[,index_test_data[[i]]],
                                     size.factor,     plot.dispersion, output.dir,
                                     remove.false )
        }

      }

    }

    res <- list()
    for (i in c(1,3,4,5,7)){
      res[[i]] <- as.data.frame(matrix(nrow = dim(meth_test)[1]))
      for (j in 1: length(index_test_data)){

        if (length(res_list[[j]])!=0){
          res[[i]] <- cbind(res[[i]],res_list[[j]][[i]])
        }

      }
      res[[i]]<-  res[[i]][, -1]
      res[[i]]<- as.matrix( res[[i]])
      res[[i]]<-  res[[i]][, order(index_test)]


    }





    for (i in c(2,6)){
      res[[i]] <- as.data.frame(matrix(nrow = dim(meth_test)[1]))
      for (j in 1: length(index_test_data)){

        if (length(res_list[[j]])!=0){
          a <- matrix(rep(res_list[[j]][[i]],(dim(res_list[[j]][[1]])[2])),nrow=
                        (dim(res_list[[j]][[1]])[1]),
                      ncol=(dim(res_list[[j]][[1]])[2]))

          res[[i]] <- cbind(res[[i]],a)
        }

      }
      res[[i]]<-  res[[i]][, -1]
      res[[i]]<- as.matrix( res[[i]])
      res[[i]]<-  res[[i]][, order(index_test)]


    }


    for (i in 1: length(res)){
      res[[i]] <- data.frame(res[[i]])
      colnames(res[[i]]) <- paste0("cell_",(1:dim(res[[i]])[2]))
    }
    res_final <-list()
    res_final [[1]] <- res[[1]]
    res_final [[2]] <- res[[2]]
    res_final [[3]] <- res[[3]]
    res_final [[4]] <- res[[4]]
    res_final [[5]] <- res[[5]]
    res_final [[6]] <- res[[6]]
    res_final [[7]] <- res[[7]]


    names(res_final )<- c("detect proportion treated","detect proportion control","log2 Risk Ratio",
                          "log2 Odds Ratio","p value","aboundance","adjusted p value")

    return(res_final )



  }
