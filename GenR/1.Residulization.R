###########################################
########### data residulization ###########
###########################################

setwd("V:/medewerkers/051950 Xu, B/PhD projects/Structural_CCA")
all_final <- readRDS("all_final_genr_zscores_30.01.2023.rds")
names(all_final)
# brain measures are 13:87
# regress the feature matrix on age, gender, ethnicity and SES (maternal education)
resi_brain <- lapply(12:120, function(i) {
  out <- residuals(lm(all_final[, i] ~ all_final$age + all_final$gender + 
                        all_final$ethnicity + all_final$maternal_edu))
})


brain_residual <- do.call(cbind, resi_brain)
#brain_residual <- cbind(IDC=all_final$IDC, brain_residual)
colnames(brain_residual) <- names(all_final)[12:120]
head(brain_residual)
saveRDS(brain_residual, "brain_residual_tbv_zscores_09.01.2023.rds")

cbcl <- all_final[,c(1:8)]
df_all <- cbind(cbcl,brain_residual)

names(df_all)
head(df_all)
saveRDS(df_all, "df_all_genr_30.01.2023.rds")