#' This function serves as the central component for single-cell RNA methylation
#' site analysis by clustering cells based on gene expression or read counts,
#' encompassing both test and control groups.
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
#' @param plot.dispersion whether plot or not
#' @param output.dir the path
#' @param remove.false remove some false positive or not
#'
#'
#' @return a list c("detect proportion treated","detect proportion control","log2 Risk Ratio",
#'  "log2 Odds Ratio","p value","abundance","adjusted p value")
#'
#'  1.The first part of the list is the meth proportion data frame (row is gene
#'  and column is single cell) in the test cells.
#'  2.The second part of the list is the mean meth proportion vector in the control cells.
#'  3.The third part of the list is the log2 risk ratio data frame for the test cells.
#'  4.The forth part of the list is the log2 odds ratio data frame for the test cells.
#'  5.The fifth part of the list is the p value data frame for the test cells.
#'  6.The sixth part of the list is the estimated gene abundance vector.
#'  7.The seventh part of the list is the adjusted p value data frame for the test cells.
#'
#' @export
#'
#' @importFrom stats p.adjust
#' @importFrom utils data write.table
#' @export
#'
#' @examples
#'   data <- simulateData (test_num=1,control_num=30)
#'
#'
#'   res <- SIGMRtest(data[[1]],data[[2]],data[[3]],data[[4]])
#'
#'   data <- simulateData (test_num=2,control_num=30)
#'
#'
#'   res <- SIGMRtest(data[[1]],data[[2]],data[[3]],data[[4]])
#'
SIGMRtest <-
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
    mean <-getEstimate(meth_control,meth_test,unmeth_control,unmeth_test,s_control,s_test)
    p0 <- mean[[1]]
    p_control <- mean[[2]]
    p_test <- mean[[3]]
    q0 <- mean[[4]]
    e_control <- mean[[5]]
    e_test <- mean[[6]]
    res <-baseFit(meth_control,unmeth_control,p_control,q0,e_control,s_control)
    fit_meth<-res[[1]]
    fit_unmeth<-res[[2]]




    if (is.na(output.dir)) {
      output.dir <- getwd()
    }

    path <- paste(output.dir,"dispersion.pdf",sep = '/')
    if(plot.dispersion){
      plotDispersion(fit_meth,fit_unmeth,path)
    }


    # calculate z
    z_t1 <-calculateZ(q0,p_control,s_control,e_control)
    z_c1 <-calculateZ(q0,(1-p_control),s_control,e_control)

    # get estimate w
    size1_t <- size2_t <- size1_c <- size2_c <- data.frame(matrix(nrow=dim(meth_test)[1],
                                                        ncol = dim(meth_test)[2]))
    mu_t1 <-mu_c1 <-mu_t2 <-mu_c2 <-data.frame(matrix(nrow=dim(meth_test)[1],ncol =
                                                        dim(meth_test)[2]))

    p0_list<- as.list(p0)



    a_function <- function (x){
      w_fit_meth <-fittedW(x,q0,fit_meth)
      w_fit_unmeth <-fittedW(x,q0,fit_unmeth)
      # get estimate of upi
      ups_t <- pmax(w_fit_meth - z_t1, 1e-8)
      ups_c <- pmax(w_fit_unmeth - z_c1, 1e-8)

      # get all variance
      raw_t1 <- raw_t2 <- rowSums(((e_control%*%t(s_control))^2*ups_t))
      raw_c1 <- raw_c2 <- rowSums(((e_control%*%t(s_control))^2*ups_c))

      mu_t1 <- rowSums(((e_control*q0*(x))%*%t(as.numeric(s_control))))
      mu_c1 <- rowSums(((e_control*q0*(1-x))%*%t(as.numeric(s_control))))


      # put size together
      size1_t <- size2_t  <- (mu_t1^2)/raw_t1
      size1_c <- size2_c  <- (mu_c1^2)/raw_c1
      list_return <- list()
      list_return [[1]] <- mu_t1
      list_return [[2]] <- mu_c1
      list_return [[3]] <- size1_t
      list_return [[4]] <- size1_c
      return(list_return)

    }
    #get estimate mean and size
    res_list <- lapply(p0_list, a_function)

    for (i in 1:dim(meth_test)[2]){
      mu_t1[,i] <- res_list[[i]][[1]]
      mu_c1[,i] <- res_list[[i]][[2]]
      size1_t[,i] <- size2_t [,i] <- res_list[[i]][[3]]
      size1_c [,i] <- size2_c [,i] <-  res_list[[i]][[4]]
    }




    mu_t2 <- t(t(e_test*q0*p0)*(s_test))
    mu_c2 <- t(t(e_test*q0*(1-p0))*(s_test))






    # observation together calculate the p value
    res <- data.frame(matrix(nrow = dim(meth_test)[1],ncol =dim(meth_test)[2]))


    for (i in 1:dim(meth_test)[2]){
      t1 <- rowSums(meth_control)
      t2 <- meth_test[,i]
      c1 <- rowSums(unmeth_control)
      c2 <- unmeth_test[,i]
      t <- t1 + t2
      n1 <- t1 + c1
      n2 <- t2 + c2
      # go to test
      res[,i] <-quadNBtest(t1,t,n2,mu_t2[,i],mu_c2[,i],
                           size2_t[,i],size2_c[,i])
      if (remove.false==TRUE){
        index <- which(p_test[,i]<=p_control)
        res[index,i] <- 1
      }
    }



    # add fc
    log2.RR <- log2(p_test/p_control)

    p.treated <- p_test
    p.control <- p_control

    log2.OR <- data.frame(matrix(nrow = dim(meth_test)[1],ncol =dim(meth_test)[2]))
    for (i in 1:dim(meth_test)[2]){
      log2.OR[,i] <- log2((rowSums(t(t(meth_test[,i])/s_test[i]))/rowSums(t(t(
        unmeth_test[,i])/s_test[i])))/(rowSums(t(t(meth_control)/s_control))/
                                         rowSums(t(t(unmeth_control)/s_control))))

    }


    padj <- as.data.frame(apply( res, 2, p.adjustfun))


    res <-list(p.treated,p.control,log2.RR,log2.OR,res,q0,padj)

    for (i in c(1,3,4,5,7)){
      res[[i]] <- data.frame(res[[i]])
      colnames(res[[i]]) <- paste0("cell_",(1:dim(res[[i]])[2]))
    }

    names(res) <-c("detect proportion treated","detect proportion control","log2 Risk Ratio",
                   "log2 Odds Ratio","p value","aboundance","adjusted p value")
    return(res)}
