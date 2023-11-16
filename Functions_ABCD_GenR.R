#######################################################################
############## important functions used in the project ################
#######################################################################


######################################################
####### Read the rs-fMRI matrices and vectorize ######
######################################################
# this is for ABCD
make_brain_features <- function(subid){
    
    conMatDir <- ('/gpfs/work2/0/einf1049/scratch/bxu/ABCD_data_filtered/')
    cl <- makePSOCKcluster(32)
    registerDoParallel(cl)
    
    feature_list <- foreach::foreach(i = 1:length(subid), .packages=c("doParallel", "foreach")) %dopar% {
        subID <- subid[i]
        conMatCsv <- file.path(conMatDir, paste(subID, 'corMat_Gordon_fsSC_filtered', sep='_'))
        conMatCsv <- paste(conMatCsv, 'csv', sep='.')
        conMat <- as.matrix(read.table(conMatCsv))
        feature_column <- conMat[upper.tri(conMat, diag = FALSE)]
    }
    
    stopCluster(cl)
    
    feature_abcd <- do.call(rbind,feature_list)
    feature_abcd <- as.data.frame(feature_abcd)
    colnames(feature_abcd) <- paste0("F_", 1:ncol(feature_abcd))
    
    return(feature_abcd)
}


#############################
####### Residualizaiton #####
#############################


residualization <- function(brain,confounders){
    cl <- makePSOCKcluster(32)
    registerDoParallel(cl)
    
    residual_list <- foreach::foreach(i = seq_along(brain), .packages=c("doParallel", "foreach")) %dopar% {
        out <- residuals(lm(brain[, i] ~ confounders$age + confounders$sex + 
                          confounders$race_ethnicity + confounders$site + confounders$parental_education, na.action=na.exclude))
    }
    
    stopCluster(cl)
    
    brain_residual <- do.call(cbind, residual_list)
    return(brain_residual)
}


###########################
####### weighted PCA ######
###########################


weighted_pca <- function(cbcl,brain,n) {
  rank <- rank(-rowSums(cbcl))
  weights <- log(nrow(cbcl)) - log(rank)
  weights <- weights / sum(weights)
  # center the original data
  feature_brain_centered <- scale(brain, scale = FALSE)
  # scale the columns of the feature matrix
  feature_centered_weighted <- feature_brain_centered * replicate(ncol(feature_brain_centered), weights)
  # call PCA on the weighted feature matrix
  pca.weighted <- prcomp(feature_centered_weighted)
  # the rotation transformation: take the first 100 eigenvectors
  rotation <- pca.weighted$rotation[,1:n]
  # the data with the reduced dimensionality
  feature_brain_reduced <- feature_brain_centered %*% rotation 
  return(brain_train_reduced=feature_brain_reduced, rotation = pca.weighted$rotation)
}


########################################
####### penalty parameters search ######
########################################


####### 1. create 100 resamples of the ABCD training set (100 training and validation sets)

CV_sampling <- function(cbcl, brain, s_num, partition) {
  
  # s_num: number of resamples
  # partition: proportion of training set
  brain_train_val <- list()
  brain_test_val <- list()
  cbcl_train_val <- list()
  cbcl_test_val <- list()
  
  for (i in 1:s_num) {
    
    df.train.validation <- data.frame(IDC=1:nrow(cbcl))
    df_train_val <- partition(data = df.train.validation, p = partition)
    
    train <- unlist(df_train_val[[1]])
    test <- unlist(df_train_val[[2]])
    
    feature_val <- brain[train, ]
    feature_val_test <- brain[test, ]
    cbcl_val <- cbcl[train, ]
    cbcl_val_test <- cbcl[test, ]
    
    if (any(colSums(cbcl_val_test) == 0)){
      brain_train_val[[i]] <- NULL
      brain_test_val[[i]] <- NULL
      cbcl_train_val[[i]] <- NULL
      cbcl_test_val[[i]] <- NULL
    }
    
    brain_train_val[[i]] <- feature_val
    brain_test_val[[i]] <- feature_val_test
    cbcl_train_val[[i]] <- cbcl_val
    cbcl_test_val[[i]] <- cbcl_val_test
    
  }
  resample <- list(brain_train=brain_train_val, brain_test=brain_test_val, 
                   cbcl_train=cbcl_train_val, cbcl_test=cbcl_test_val)
  return(resample)
}



######## 2. Penalty parameters selection:part of the code is based on Xia et al.(2018) 

cca_resample_test_m <- function(X, Y, X2, Y2, pen_x, pen_y,nsample) { 
  # X and Y are training sets created by CV_sampling, 
  # X2 and Y2 are test sets created by CV_sampling
  # pen_x is for brain
  # pen_y is for cbcl
  # nsample is how many resamples you have 
  
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  cca.out <- foreach::foreach(i=1:nsample, .packages=c("PMA","doParallel","permute")) %dopar% {
    cca.res <- function(X, Y, pen_x, pen_y, cv_num) {
      cca.out <- PMA::CCA(x=X, z=Y, typex="standard", typez="standard", penaltyx = pen_x, penaltyz = pen_y, 
                          niter = 20, K=cv_num)
    }
    # calculate the correlations in the test sets
    res <- cca.res(X[[i]], Y[[i]], pen_x, pen_y, 8)
    std_brain <- scale(X2[[i]]) %*% res$u
    std_cbcl <- scale(Y2[[i]]) %*% res$v
    
    cor.res <- diag(cor(std_brain, std_cbcl))
  }
  stopCluster(cl)
  # extract the first canonical correlations
  cca.cor <- sapply(1:length(X), function(i) cca.out[[i]][1])
  cor.mean <- mean(abs(cca.cor), na.rm=TRUE)
  list(Penalty_x=pen_x, Penalty_y=pen_y,Cor_mean=cor.mean)
}


######## 3. grid search: code is based on Xia et al. (2018)
                    
grid.search.cor.Testset <- function(X,Y,X2,Y2,pen_xseq,pen_yseq,nsample) {
  # X and Y are training sets created by CV_sampling, 
  # X2 and Y2 are test sets created by CV_sampling
  # pen_xseq is for brain
  # pen_yseq is for cbcl
  # nsample is how many resamples you have
  
  # loop through x penalty
  pen_x_loop <- function(X,Y,X2,Y2,pen_xseq,pen_y,nsample){
    pen_x_loop <- lapply(pen_xseq, function(x) cca_resample_test_m(X,Y,X2,Y2,x,pen_y,nsample))
  }
  # then loop through y penalty
  pen_y_loop <- lapply(pen_yseq, function(y) pen_x_loop(X,Y,X2,Y2,pen_xseq,y,nsample))
  
  cor.xy <- unlist(pen_y_loop, recursive = F)
  cor.res <- t(sapply(cor.xy, function(x) unlist(x)))
  best.para <- cor.res[which.max(abs(cor.res[, "Cor_mean"])), ]
  
  cor.mat <- matrix(cor.res[, "Cor_mean"], 10, 10)
  rownames(cor.mat) <- x_pen
  colnames(cor.mat) <- y_pen
  return(list(best.para, cor.mat))
}  


###################################
####### covariance explained ######
###################################

###### Code is based on Xia et al.(2018)
VarianceExplain <- function(brain,cbcl,ccares, n) {
  residual_std <- apply(brain, 2, scale)
  cbcl_std <- apply(cbcl, 2, scale)
  
  covmat <- t(ccares$u) %*% t(residual_std) %*% cbcl_std %*% ccares$v
  varE <- diag(covmat)^2 / sum(diag(covmat)^2) 
  varE.df <- data.frame(1:n, var=varE)
  
  plot(varE.df)
  abline(h=mean(varE.df$var), col="red", lty=3, lwd=2)
  return(varE.df)
}

                     
                      
####################################
###### Generalizability test  ######
####################################


# project the weights and calculate canonical correlations
test_project_weights <- function(brain,cbcl,training_model, nCor){
  std_brain <- scale(brain) %*% training_model$u
  std_cbcl <- scale(cbcl) %*% training_model$v
  cor <- round(diag(cor(std_brain, std_cbcl)), nCor)
  return(cor)
}


# permutation tests in ABCD training sets
permutation_test <- function(cbcl,brain,nperm, penaltyCBCL,penaltyBrain,numCor, cors) { 
  # numCor is the number of correlations you are interested in
  # cors are the correlations in the unshuffled data
  n.perm = nperm
  shuffle_idx <- sapply(1:n.perm, function (x){permute::shuffle(1:nrow(cbcl))})
  cbcl_perm <- lapply(1:n.perm, function(i) {cbcl[shuffle_idx[, i], ]})
  
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  perm_cor <- foreach::foreach(i = seq_along(cbcl_perm)) %dopar% {
    library(PMA)
    res_perm <- CCA(x=cbcl_perm[[i]], z=brain, penaltyx = penaltyCBCL,  
                    penaltyz = penaltyBrain, niter=20, K=numCor)
  }
  stopCluster(cl)
  
  cor.perm <- lapply(perm_cor, function(x) {x$cors})
  cor.perm <- do.call(rbind, cor.perm)
  pval.perm <- sapply(1:numCor, function(x){length(which(abs(cor.perm[, x]) >= abs(cors[x])))/(n.perm+1)})
  
  out <- list(cor.perm=cor.perm, pval.perm=pval.perm)
  return(out)
  
}

# permutation tests in ABCD test sets
permutation_test_testset <- function(cbcl,brain,nperm, model,cors) { 
  n.perm = nperm
  shuffle_idx <- sapply(1:n.perm, function (x){permute::shuffle(1:nrow(cbcl))})
  cbcl_perm <- lapply(1:n.perm, function(i) {cbcl[shuffle_idx[, i], ]})
  
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  perm_cor <- foreach::foreach(i = seq_along(cbcl_perm)) %dopar% {
    std_brain <- scale(brain) %*% model$u
    std_cbcl <- scale(cbcl_perm[[i]]) %*% model$v
    cor <- round(diag(cor(std_brain, std_cbcl)), 6)
  }
  stopCluster(cl)
  
  cor.perm <- do.call(rbind, perm_cor)
  pval.perm <- sapply(1:6, function(x){length(which(abs(cor.perm[, x]) >= abs(cors[x])))/(n.perm+1)})
  
  return(list(cor.perm=cor.perm, pval.perm=pval.perm))
}













