######################################################################################
########### GenR data preparation: brain structures, covariates, and CBCL ############
######################################################################################

library(data.table)
library(dplyr)
library(stringr)
library(Amelia)

###########################
#structural brain data QC
###########################
setwd("V:/medewerkers/***** Xu, B/PhD projects/Structural_CCA")

# read the QC data of brain measures
core <- readRDS("genr_mri_core_data_20220311.rds")

###### 1. filter the ids that passed the brain imaging QC
qc <- core %>% filter(mri_consent_f09 == "yes" & t1_has_nii_f09 == "yes" & 
                      t1_asset_has_nii_f09 != "exclude" & has_braces_mri_f09 == "no" &
                      exclude_incidental_f09 == "include" & freesurfer_qc_f09 == "usable")

# then selected the covariate (age at the MRI) that we are interested in. 
core_qc <- qc[, c("idc","age_child_mri_f09")]
dim(core_qc)

###### 2. remove twins/siblings
# read the twin/siblings information of children.
general <- read.csv("CHILD-ALLGENERALDATA_12112020.csv")
# change the idc name for data merging.
names(general)[1] <- "idc"

###### 3. merge the data
qc_all <- merge(core_qc, general[,c("idc","MOTHER")], by="idc")
dim(qc_all)
# randomly included one twin/siblings
qc_twins <- qc_all[!duplicated(qc_all$MOTHER),]
dim(qc_twins)

###########################
#structural data extraction
###########################
# this process is the same with ABCD

cortical <- readRDS("f09_freesurfer_v6_09dec2016_aparc_stats_pull06june2017.rds")
cortical <- cortical[cortical$idc %in% qc_twins$idc, ]
subcortical <- readRDS("f09_freesurfer_v6_09dec2016_aseg_stats_pull06june2017_v1.rds")
subcortical <- subcortical[subcortical$idc %in% qc_twins$idc, ]
names(cortical)
names(subcortical)

# function of average left and right cortical regions
ave_lhrh <- function(measures) {
  cor <- cortical %>% select(contains(measures)) %>% as.data.frame
  names(cor) <- names(cor) %>% gsub(c("lh_|rh_|_f09"),"",.)
  # loop through all the cortical regions, there are 34 regions
  # and then average left and right
  temp1 <- lapply(1:34, function(x) rowMeans(cor[,names(cor)==names(cor)[x]])) # the region names
  temp2 <- do.call(cbind, temp1)
  colnames(temp2) <- names(cor)[1:34]
  temp2 <- as.data.frame(temp2)
  return(temp2)
}


######## cortical volumes
cor_vol <- ave_lhrh("vol")

######## cortical thickness
cor_thick <- ave_lhrh("thickavg")

######## cortical reas 
cor_surfarea <- ave_lhrh("surfarea")

######## total brain volumes
tbv <- global[,"genr_tbv_f09"]


######## subcortical volumes

# subcortical areas of interest
sub_area <- c("Hippocampus","Amygdala","Accumbens","Caudate","Thalamus","Putamen","Pallidum")

sub1 <- lapply(sub_area, function(x) rowMeans(subcortical[, str_detect(names(subcortical),x)]))# the region names
sub2 <- do.call(cbind, sub1)
colnames(sub2) <- paste0(sub_area, "_vol")
subcor_vol <- as.data.frame(sub2)

###########################
#CBCL syndrome scores extraction
##########################

cbcl_total <- fread("CHILDCBCL9_incl_Tscores_20201111.csv")
names(cbcl_total)

cbcl_syndrom <- cbcl_total[cbcl_total$IDC %in% qc_twins$idc, 
                           c("IDC","sum_anx_9m","sum_wit_9m","sum_som_9m","sum_sop_9m",
                               "sum_tho_9m","sum_att_9m","sum_rul_9m","sum_agg_9m", 
                               "sum_int_9m","sum_ext_9m","sum_oth_9m","cbcl_sum_9m")]

cbcl_syndrom[cbcl_syndrom == 888] <- NA
cbcl_syndrom[cbcl_syndrom == 999] <- NA

# combine all the data
cbcl_brain <- cbind(cbcl_syndrom, cor_vol, cor_thick, cor_surfarea, subcor_vol)
str(cbcl_brain)
dim(cbcl_brain)


###########################
#covariates extraction
###########################
# read the covariates we are interested: parental nation of origin, sex, maternal education, family income was used for imputation. 
demo <- general[,c("idc","IDM","ETHNINFv3","GENDER","EDUCM5","INCOME5")]
               
# merge all the data
df_all <- merge(cbcl_brain, demo, by="idc")
df_all <- merge(df_all, qc_twins, by="idc")
dim(df_all)
names(df_all)
summary(df_all)

###########################
#deal with NAs
###########################

# check all the NAs
colSums(is.na(df_all))

# remove the children with > 25% missingness in CBCL
# there are 8 syndrome scores, so > 8 * 75% = 6
df_all[df_all == 999] <- NA
df_all[df_all == 888] <- NA
df_all[df_all == 777] <- NA

df_all_nacbcl <- df_all[rowSums(is.na(df_all[,c("sum_anx_9m","sum_wit_9m","sum_som_9m","sum_sop_9m",
                               "sum_tho_9m","sum_att_9m","sum_rul_9m","sum_agg_9m")])) < 6, ]

# check all the NAs
colSums(is.na(df_all_nacbcl))

###########################
#EM imputation
###########################

df_imp <- df_all_nacbcl %>% mutate(age = scale(age_child_mri_f09), gender=GENDER, 
                                   ethnicity=factor(ETHNINFv3), income=factor(INCOME5),
                                   maternal_edu=factor(EDUCM5))
# remove info that are not useful
dfimp1 <- df_imp[, !names(df_imp) %in% c("cbcl_sum_9m","IDC","IDM","ETHNINFv3",
                           "GENDER","EDUCM5","INCOME5","age_child_mri_f09")]
str(dfimp1)
               
# EM imputation
all_toimp <- amelia(dfimp1, noms=c("ethnicity"),ords=c("maternal_edu","income"),
                    m=1, boot.type = "none")

all_imp <- cbind(idc=df_all_nacbcl$idc,all_toimp$imputations$imp1)

# change the levels of ethnicity & maternal education
levels(all_imp$ethnicity) <- list("dutch" = 1,
                                  "european" = c(7,9,700),
                                  "others" = c(2:6,8,10,200,300,400,500,600,800))

levels(all_imp$maternal_edu) <- list("low" = 0:1, 
                                     "medium" = 2:3,
                                     "high" = 4:5)


###########################
#standardize the brain measures
###########################
# final step: standardize the brian measures

# select the brain measures: cortical volumes, surface areas, thickness, subcortical volumes
all_imp[13:121] <- lapply(all_imp[13:121], scale)
# remove ids
all_genr_zscores <- all_imp[,-1]
dim(all_genr_zscores)

saveRDS(all_genr_zscores , "all_final_genr_zscores.rds")
