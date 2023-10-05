########################################################
####################### ABCD ##########################
########################################################


########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'fmsb', 'NbClust','stats','dendextend','cluster',
              'fpc', 'e1071', 'plotly', 'MASS', 'ggradar2','tidyr','introdatavizyes','Hmisc',
              'pheatmap','RColorBrewer','circlize','CCA','CCP','candisc')


lapply(packages, require, character.only = TRUE)


########################################################
############## 1. Read the data ########################
########################################################
all_abcd <- readRDS("train_test_split_abcd.rds") # this is the data from previous step: Residualization_abcd.R

########################################################
############## 2. split in training and test ###########
########################################################
### read the training data and test data

str_train <- lapply(1:10, function(i) {temp <- all_abcd[[i]]
                                      list(brain_train=temp$brain_train,
                                           cbcl_train=temp$cbcl_train)}) 

str_test <- lapply(1:10, function(i) {temp <- all_abcd[[i]]
                                      list(brain_test=temp$brain_test,
                                           cbcl_test=temp$cbcl_test)})

########################################################
############## 3. Grid search for penalty parameters  ##
########################################################
grid.abcd.str.all <- lapply(1:10, function(i) {
  structure_resample <- str_train[[i]]
  brain_train <- structure_resample$brain_train
  cbcl_train <- structure_resample$cbcl_train
  cv.resample <- CV_sampling(cbcl_train, brain_train, 100, 0.8)
  
  x_pen <- seq(0.1, 1.0, length.out = 10) # brain
  y_pen <- seq(0.1, 1.0, length.out = 10)

  grid.abcd <- grid.search.cor.Testset(cv.resample$brain_train, cv.resample$cbcl_train,
                                        cv.resample$brain_test, cv.resample$cbcl_test,
                                        x_pen, y_pen, nsample=100)
  })

saveRDS(grid.abcd.str.all,"grid.abcd.str.all.rds")


###################################################################################################
############## 4. Fit the sCCA model  across ABCD training, test and Genertion R ##################
###################################################################################################

str_train_test_abcd <- lapply(1:10, function(i) {

  ######### first extract the data from each train-test split
  
  brain_train <- str_train[[i]]$brain_train
  cbcl_train <- str_train[[i]]$cbcl_train
  brain_test <- str_test[[i]]$brain_test
  cbcl_test <- str_test[[i]]$cbcl_test
  brain_pen <- grid.abcd.str.all[[i]][[1]][1]
  cbcl_pen <- grid.abcd.str.all[[i]][[1]][2]

  ######### fit the SCCA model in ABCD training set and do permuation test
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  perm_abcd_train <- permutation_test(cbcl_train, brain_train, nperm=999, cbcl_pen, brain_pen, 8, res.abcd$cors)

  ######### project the weights from ABCD training set to ABCD test set and do permutation test
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd, 8)
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=1999,res.abcd,abcd.test)

  ######### project the weights from ABCD training set to Generation R and do permutation test
  cor.abcdTogenr <- test_project_weights(brain_genr_str, cbcl_genr_str, res.abcd,8)
  perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr_str,brain_genr_str,nperm=1999, 
                                                    res.abcd,abs(cor.abcdTogenr))

  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr_total$pval.perm))
})


saveRDS(str_train_test_abcd,"str_train_test_abcd_genr.rds")





