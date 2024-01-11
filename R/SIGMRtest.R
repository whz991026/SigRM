#' The main test function
#' 
#' 
#' 
#' @description This function serves as the central component for single-cell RNA methylation
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
#' @param q expression aboundance
#' @param method_dispersion the method to estimate the dispersion ("unbiased_estimate","mle","locfit")
#' @param one_side p value is one side or both sides
#' @param adjust_1 adjust parameter for RR
#' @param adjust_2 adjust parameter for OR
#' 
#' @return a list c("methylation level treated","methylation level control","log2 Risk Ratio",
#' "log2 Odds Ratio","p value","aboundance","adjusted p value","TCR")
#'
#'  1.The first part of the list is the meth proportion data frame (row is gene
#'  and column is single cell) in the test cells.
#'  2.The second part of the list is the mean meth proportion vector in the control cells.
#'  3.The third part of the list is the log2 risk ratio data frame for the test cells.
#'  4.The forth part of the list is the log2 odds ratio data frame for the test cells.
#'  5.The fifth part of the list is the p value data frame for the test cells.
#'  6.The sixth part of the list is the estimated gene abundance vector.
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
           remove.false = TRUE,
           q=NA,
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
    
    if(  any(is.na(s_control)) | any(is.na(s_test))){
      stop( "size factor contains NA" )
    }
    # estimate
    print("Estimating dispersion for each RNA methylation site, this will take a while ...")
    
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
    
 
    
    if(method_dispersion=="locfit"){
      res <-baseFit(meth_control,unmeth_control,p_control,q0,e_control,s_control)
      fit_meth<-res[[1]]
      fit_unmeth<-res[[2]]
     
      a <- try(fitted_v(0.5,100,fit_meth))
      b <- try(fitted_v(0.5,100,fit_unmeth))
      
      
      if(inherits(a,"try-error")){
        fit_meth <- FALSE
        print("locfit for meth reads get error replace to unbiased estimate")
      }
      
      if(inherits(b,"try-error")){
        fit_unmeth <- FALSE
        print("locfit for unmeth reads get error replace to unbiased estimate")
      }
   
        if (is.na(output.dir)) {
          output.dir <- getwd()
        }
        
        path <- paste(output.dir,"dispersion.pdf",sep = '/')
        if(plot.dispersion){
          if(inherits(fit_meth,"locfit")|inherits(fit_meth,"locfit")){
            plotDispersion(fit_meth,fit_unmeth,path)
          }else{
            print("locfit for both get error can not plot")
          }
          
        }
      
      
    }
    
    
    
    if(method_dispersion=="locfit"|method_dispersion=="unbiased_estimate"){
      
      
      # get estimate v
      
      
      p0_list<- as.list(p0)
      
      
      
      a_function <- function (x){
        
        if(method_dispersion=="unbiased_estimate"){
          v_fit_U <-as.numeric(calculate_v (meth_control,s_control,e_control))
          v_fit_C <-as.numeric(calculate_v (unmeth_control,s_control,e_control))
        }else{
          
          if(inherits(fit_meth,"locfit")){
            v_fit_U <-fitted_v(x,q0,fit_meth)
          } else{
            v_fit_U <-as.numeric(calculate_v (meth_control,s_control,e_control))
          }
          
          if(inherits(fit_unmeth,"locfit")){
            v_fit_C <-fitted_v((1-x),q0,fit_unmeth)
          } else{
            v_fit_C <-as.numeric(calculate_v (unmeth_control,s_control,e_control))
          }
          
        }
        
        
        # calculate z
        z_U_control <-calculateZ(q0,x,s_control,e_control)
        z_C_control <-calculateZ(q0,(1-x),s_control,e_control)
        
        # get estimate of upi
        ups_U <- pmax(v_fit_U - z_U_control, 1e-8)
        ups_C <- pmax(v_fit_C - z_C_control, 1e-8)
        
        # get all variance
        raw_U <- (((e_control%*%t(s_control))^2)*ups_U)
        raw_C <- ((e_control%*%t(s_control))^2*ups_C)
        
        mu_U <- ((e_control*q0*(x))%*%t(as.numeric(s_control)))
        mu_C <- ((e_control*q0*(1-x))%*%t(as.numeric(s_control)))
        
        
        # put size together
        size_U <- ((mu_U^2)/raw_U)[,1]
        size_C <- ((mu_C^2)/raw_C)[,1]
        list_return <- list()
        list_return [[1]] <- size_U
        list_return [[2]] <- size_C
        return(list_return)
        
      }
      #get estimate mean and size
      res_list <- lapply(p0_list, a_function)
      
      size_U <- as.data.frame(lapply(res_list,function(x) x[][[1]]))
      size_C <- as.data.frame(lapply(res_list,function(x) x[][[2]]))
      
      
      
    }else{
      mu_U_control <- as.matrix(e_control*q0*p_control)%*%(s_control)
      mu_C_control <- as.matrix(e_control*q0*(1-p_control))%*%(s_control)
      
      
      
      
      size <- Estimate_dispersion(meth_control,unmeth_control,mu_U_control,mu_C_control)
      
      
      index_na <- which(is.na(mu_U_control[,1]))
      size_U <- as.numeric(size[[1]])
      size_C <- as.numeric(size[[2]])
      
      size_U[index_na] <- NA
      size_C[index_na] <- NA
      
      
      
    }
    
    



    mu_U_test_0 <- t(t(e_test*q0*p0)*(s_test))
    mu_C_test_0 <- t(t(e_test*q0*(1-p0))*(s_test))


    





    # observation together calculate the p value
    
    if(method_dispersion=="mle"){
      matrix_num <- matrix(rep(dim(meth_test)[1], dim(meth_test)[2]),nrow = 1)
      colnames(matrix_num) <- colnames(unmeth_test)
      
      
      list_create <- as.list(as.data.frame(rbind(matrix_num,as.matrix(meth_test), 
    as.matrix(unmeth_test),as.matrix(mu_U_test_0),as.matrix(mu_C_test_0),
    matrix(rep(size_U,dim(meth_test)[2]),
   ncol=dim(meth_test)[2]),matrix(rep(size_C,dim(meth_test)[2]),
      ncol=dim(meth_test)[2]))))
      
    }else{
      
      matrix_num <- matrix(rep(dim(meth_test)[1], dim(meth_test)[2]),nrow = 1)
      colnames(matrix_num) <- colnames(unmeth_test)
      
      
      list_create <- as.list(as.data.frame(rbind(matrix_num,as.matrix(meth_test),
                     as.matrix(unmeth_test),as.matrix(mu_U_test_0),
                     as.matrix(mu_C_test_0),as.matrix(size_U),as.matrix(size_C))))
    }
    
    


    list_function <- function(x){
      x <- unlist(x)
      t2 <- x[2:(1+x[1])]  
      c2 <- x[(2+x[1]):(1+2*x[1])] 
      n2 <- t2 + c2
      
      
      
      res <-quadNBtest(t2,n2,x[(2+2*x[1]):(1+3*x[1])],x[(2+3*x[1]):(1+4*x[1])],
                       x[(2+4*x[1]):(1+5*x[1])],x[(2+5*x[1]):(1+6*x[1])],one_side=one_side)
      
      return(res)
    }
    
    
    res <- as.data.frame(lapply(list_create, list_function))
    
    
   
    
    
    
    p.treated <- p_test
    p.control <- p_control
    
    # add RR 
    log2.RR <- log2((p_test+adjust_1)/(p_control+adjust_1))
    
    
    if (remove.false==TRUE){
      res <- as.list(res)
      list_create <- as.list(1:dim(meth_test)[2])
      list_function <- function(x){
        index <- which(log2.RR[,x]<0)
        res[[x]][index] <- 1
        res[[x]]
      }
      res <- as.data.frame(lapply(list_create, list_function))
    }
    
    
    # add OR 
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
    
    
    
    padj <- as.data.frame(apply( res, 2, p.adjustfun))
    
    TCR <- p.treated-p.control
    
    TCR <- data.frame(TCR)
    
    p.treated <- data.frame(p.treated)
    colnames(p.treated) <- paste0("cell_",(1:dim(p.treated)[2]))
    
    log2.RR <- data.frame(log2.RR)
    colnames(log2.RR) <- paste0("cell_",(1:dim(log2.RR)[2]))
    
    log2.OR <- data.frame(log2.OR)
    colnames(log2.OR) <- paste0("cell_",(1:dim(log2.OR)[2]))
    
    res <- data.frame(res)
    colnames(res) <- paste0("cell_",(1:dim(res)[2]))
    
    padj <- data.frame(padj)
    colnames(padj) <- paste0("cell_",(1:dim(padj)[2]))
    
    TCR <- data.frame(TCR)
    colnames(TCR) <- paste0("cell_",(1:dim(TCR)[2]))
    
    res <-list(p.treated,p.control,log2.RR,log2.OR,res,q0,padj,TCR)
    
    
    
    names(res) <-c("methylation level treated","methylation level control","log2 Risk Ratio",
                   "log2 Odds Ratio","p value","aboundance","adjusted p value","TCR")
return(res)}