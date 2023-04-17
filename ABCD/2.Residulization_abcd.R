###########################################
########### data residulization ###########
###########################################

setwd("/home/r051950/structural_cca/abcd/")
all_final <- readRDS("all_final_abcd.rds")

residualization <- function(brain,confounders){

     resi_brain <- lapply(1:ncol(brain), function(i) {
        out <- residuals(lm(brain[, i] ~ confounders$interview_age + confounders$sex + 
        confounders$race_ethnicity + confounders$site + confounders$parental_education, na.action=na.exclude))
     })
     
     brain_residual <- do.call(cbind, resi_brain)
     return(brain_residual)
}
    
############### read the data
train <- readRDS("all_final_train_tbv.rds")
test <- readRDS("all_final_test_tbv.rds") 

train_test_split <- lapply(1:10, function(i) {

    train0 <- train[[i]]
    subid_train <- train0$idc
    subid_train <- as.character(subid_train)
    
    test0 <- test[[i]]
    subid_test <- test0$idc
    subid_test <- as.character(subid_test)
    
    ########### data extraction
    brain_train <- all_final[all_final$idc %in% subid_train, 23:131] 
    brain_train <- apply(brain_train, 2, scale)
    brain_test <- all_final[all_final$idc %in% subid_test, 23:131]
    brain_test <- apply(brain_test, 2, scale)
    
    #cbcl_train <- train0[,10:128]
    #cbcl_test <- test0[,10:128]
    
    cbcl_train <- train0[, c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                          "cbcl_scr_syn_social_r", "cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                          "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")]
    cbcl_test <- test0[,c("cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                          "cbcl_scr_syn_social_r", "cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                          "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r")]
                           
    #confounders_train <- train0[,c("interview_age","sex","race_ethnicity","site","parental_education","SurfaceHoles")]
    #confounders_test <- test0[,c("interview_age","sex","race_ethnicity","site","parental_education","SurfaceHoles")]
    
    confounders_train <- all_final[all_final$idc %in% subid_train,c("interview_age","sex","race_ethnicity","site","parental_education","SurfaceHoles")]
    confounders_test <- all_final[all_final$idc %in% subid_test,c("interview_age","sex","race_ethnicity","site","parental_education","SurfaceHoles")]
    
    ############ residualization
    brain_train_residual <- residualization(brain_train, confounders_train)
    colnames(brain_train_residual) <- names(all_final)[23:131]
    brain_test_residual <- residualization(brain_test, confounders_test)
    colnames(brain_test_residual) <- names(all_final)[23:131]
    
    out <- list(brain_train=brain_train_residual, cbcl_train=cbcl_train, cbcl_test=cbcl_test, brain_test=brain_test_residual)

})

saveRDS(train_test_split, "train_test_split_abcd.rds")
