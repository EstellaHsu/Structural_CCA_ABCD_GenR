#########################################
############## ABCD QC ##################
#########################################


setwd("/Users/estella/Desktop/ABCD_download")

#########################################
#demographic info
#########################################
participants <- read.delim(file = 'abcd_lpds01.txt')
names(participants)
head(participants)

# rename idc for merging data
participants$idc <- sapply(as.character(participants$src_subject_id), 
                           function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

# only 1 year follow up had the parental education data
prnt_ed <- participants[-1,c("idc","demo_prnt_ed_v2_l","eventname")]
prnt_ed1 <- prnt_ed %>% filter(as.character(eventname) == "1_year_follow_up_y_arm_1")
#table(prnt_ed$demo_prnt_ed_v2_l[prnt_ed$eventname == "1_year_follow_up_y_arm_1"])
table(prnt_ed$eventname)

# other data we need baseline
participants1 <- participants %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
participants2 <- participants1[, c("idc","sex","interview_age")]

# merge the useful variables
demo_info <- merge(prnt_ed1[,c("idc","demo_prnt_ed_v2_l")],participants2,by="idc")

#########################################
#Site info
#########################################

site <- read.delim(file = 'abcd_lt01.txt')
site$idc <- sapply(as.character(site$src_subject_id), 
                   function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
site1 <- site %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
site_info <- site1[, c("idc","site_id_l")]


#########################################
#genetic: twins & siblings
#########################################

twin <- read.delim("acspsw03.txt")
names(twin)

twin$idc <- sapply(as.character(twin$src_subject_id), function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
twin1 <- twin %>% filter(eventname == "baseline_year_1_arm_1")
twin_info <- twin1[,c("idc","race_ethnicity","rel_family_id")]


#########################################
#brain data QC
#########################################

setwd("/Users/estella/Desktop/ABCD_download/manual_qc")
qc_str <- read.delim("abcd_imgincl01.txt")
names(qc_str)

qc_str$idc <- sapply(as.character(qc_str$src_subject_id), 
                               function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
qc1 <- qc_str %>% filter(eventname == "baseline_year_1_arm_1")

strqc_info1 <- qc1[,c("idc","imgincl_t1w_include")]


#########################################
#brain data incidental findings
#########################################
setwd("/Users/estella/Desktop/ABCD_download")
qc_findings <- read.delim("abcd_mrfindings02.txt")
names(qc_findings)

qc_findings$idc <- sapply(as.character(qc_findings$src_subject_id), 
                     function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
qc2 <- qc_findings %>% filter(eventname == "baseline_year_1_arm_1")

strqc_info2 <- qc2[,c("idc","mrif_score")]
str(strqc_info2)

#########################################
#merge all the data
#########################################

#put all data frames into list
df_list <- list(strqc_info1, strqc_info2, twin_info, site_info, demo_info)

#merge all data frames in list
all_info <- df_list %>% reduce(full_join, by='idc')
dim(all_info)

# 1. id has str data
setwd("/Users/estella/Desktop/Structural_cca")
id <- readRDS("id_has_strdata.rds")
all_info <- all_info %>% filter(idc %in% id)
dim(all_info)

# 1. inclusion criteria
all_info1 <- all_info %>% filter(imgincl_t1w_include == 1 & mrif_score < 3)
dim(all_info1)

# 2. remove twin/siblings
all_info2 <- all_info1[!duplicated(all_info1$rel_family_id), ]
dim(all_info2)


#########################################
#merge CBCL
#########################################

setwd("/Users/estella/Desktop/ABCD_download/ABCD_CBCL")
cbcl <- read.delim("abcd_cbcls01.txt")
cbcl$idc <- sapply(as.character(cbcl$src_subject_id),function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

cbcl1 <- cbcl %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
names(cbcl1)
######## extract syndrome scale 
cbcl_syn <- cbcl1[, c("idc","cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                      "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                      "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                      "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r")]
dim(cbcl_syn)


###### merge all the data

all_cbcl <- merge(all_info2, cbcl_syn,by="idc")
str(all_cbcl)

all_cbcl1 <- all_cbcl %>% mutate_at(c("interview_age","cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                         "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                         "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                         "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r"),as.numeric) 

colSums(is.na(all_cbcl1))

# there are many NAs in the final data in parental education, sex, interview age
# so try to read the ABCC participants data
setwd("/Users/estella/Desktop/ABCD_download")
dd_rsmri <- read.table(file = 'participants.tsv', sep = '\t', header = TRUE)
names(dd_rsmri)[1] <- "idc"
# use the parental edu, sex, and age from this dataset. 
rsmri <- dd_rsmri[, c("idc","sex","age","parental_education","site")]

apply(rsmri,2,table)
# fill the missingness
na_demo <- all_cbcl1[is.na(all_cbcl1$demo_prnt_ed_v2_l), ]
dim(na_demo)
na_fill <- merge(rsmri, na_demo[,!names(na_demo) %in% c("sex","interview_age","demo_prnt_ed_v2_l")], by=c("idc","site"))
dim(na_fill)
sum(duplicated(na_fill$idc))

all_final0 <- all_cbcl1 %>% 
  mutate(demo_prnt_ed_v2_l=replace(demo_prnt_ed_v2_l, is.na(demo_prnt_ed_v2_l), na_fill$parental_education)) %>%
  mutate(interview_age=replace(interview_age, is.na(interview_age), na_fill$age)) %>%
  mutate(sex=replace(sex, is.na(sex), na_fill$sex))

dim(all_final0)
all_final0$sex[all_final0$sex == 1] <- "M"
all_final0$sex[all_final0$sex == 2] <- "F"

#cc <- all_final0[,c("idc","interview_age","sex","demo_prnt_ed_v2_l")]
#cbind(cc[40:70,], all_cbcl1[40:70,c("idc","interview_age","sex","demo_prnt_ed_v2_l")])
#na_fill[na_fill$idc == "sub-NDARINV05WK8AN7", ]

names(all_final0)[6:7] <- c("site","parental_education")

all_final0[all_final0 == 777] <- NA
all_final0[all_final0 == 888] <- NA
all_final0[all_final0 == 999] <- NA
all_final0[all_final0 == ""] <- NA


apply(all_final0,2,table)
colSums(is.na(all_final0))
sum(duplicated(all_final0$idc))

# remove data with NAS in relevant variables
idx_noNA <- which(rowSums(is.na(all_final0[, !names(all_final0) %in% c("imgincl_t1w_include","mrif_score","rel_family_id")])) == 0)

all_final_noNA <- all_final0[idx_noNA, ]
dim(all_final_noNA)

# change the data type
str(all_final_noNA)

all_final_noNA <- all_final_noNA %>% mutate_at(c("sex","site","race_ethnicity","parental_education"), as.factor)
str(all_final_noNA)
colSums(is.na(all_final_noNA))

# adjust the level of parental education
levels(all_final_noNA$parental_education) <- list("low" = 1:12,"medium" = 13:17,"high" = 18:21)
str(all_final_noNA)
dim(all_final_noNA)

saveRDS(all_final_noNA,"all_final_noNA_update.rds")


########################################
#Add ICV, total brain volumes
#########################################
setwd("/Users/estella/Desktop/Structural_cca/")

icv <- read.delim("abcd_smrip10201.txt")
names(icv)
icv$idc <- sapply(as.character(icv$src_subject_id), 
                          function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))
icv2 <- icv %>% filter(eventname == "baseline_year_1_arm_1")

icv_tbv <- icv2[,c("idc","smri_vol_scs_intracranialv","smri_vol_scs_suprateialv")]

str(icv_tbv)


all_final_noNA <- merge(all_final_noNA, icv_tbv, by="idc")
dim(all_final_noNA)
str(all_final_noNA)

all_final_noNA <- all_final_noNA %>% mutate_at(c("smri_vol_scs_intracranialv","smri_vol_scs_suprateialv"), as.numeric)
str(all_final_noNA)
colSums(is.na(all_final_noNA))

saveRDS(all_final_noNA,"all_final_abcd_tbv.rds")

########################################
#train test split: site
#########################################

set.seed(519)
# first split the site 10 times
sites <- unique(all_final_noNA$site)
# there are 22 sites, so 22*0.8=17.6 for the training sites
train_sites <- lapply(1:10, function(i) sample(sites, size = 17, replace = FALSE))

# then split of the ids
all_final_train <- lapply(1:10, function(i) all_final_noNA %>% filter(site %in% train_sites[[i]]))
all_final_test <- lapply(1:10, function(i) all_final_noNA %>% filter(!site %in% train_sites[[i]]))

lapply(all_final_train, nrow)
names(all_final_test[[1]])
class(all_final_train[[1]])


saveRDS(all_final_train,"all_final_train_tbv.rds")
saveRDS(all_final_test,"all_final_test_tbv.rds")

setwd("/Users/estella/Desktop/Structural_cca/")

all_final_train <- readRDS("all_final_train_tbv.rds")
all_final_test <- readRDS("all_final_test_tbv.rds")

train <- all_final_train[[1]]
test <- all_final_test[[1]]

summary(test)
apply(test,2,sd)
##############################################################################################################################3
#####################################################################################
################# merge CBCL items
#####################################################################################

setwd("/Users/estella/Desktop/ABCD_download/ABCD_CBCL")
cbcl <- read.delim("abcd_cbcl01.txt")
cbcl$idc <- sapply(as.character(cbcl$src_subject_id),function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

cbcl1 <- cbcl %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
names(cbcl1)
######## extract syndrome scale 
cbcl_item <- cbcl1[, c(132,10:128)]
dim(cbcl_item)


###### merge all the data

all_cbcl_item <- merge(all_info2, cbcl_item,by="idc")
str(all_cbcl_item)
names(all_cbcl_item)

all_cbcl_item[,10:128] <- apply(all_cbcl_item[,10:128],2,as.numeric)

colSums(is.na(all_cbcl_item))

# there are many NAs in the final data in parental education, sex, interview age
# so try to read the ABCC participants data
setwd("/Users/estella/Desktop/ABCD_download")
dd_rsmri <- read.table(file = 'participants.tsv', sep = '\t', header = TRUE)
names(dd_rsmri)[1] <- "idc"
# use the parental edu, sex, and age from this dataset. 
rsmri <- dd_rsmri[, c("idc","sex","age","parental_education","site")]
names(all_cbcl_item)[6] <- "site"

all_final <- merge(rsmri, all_cbcl_item[,!names(all_cbcl_item) %in% c("sex","interview_age","demo_prnt_ed_v2_l")], by=c("idc","site"))

all_final[all_final == 777] <- NA
all_final[all_final == 888] <- NA
all_final[all_final == 999] <- NA
all_final[all_final == ""] <- NA

colSums(is.na(all_final))
sum(duplicated(all_final$idc))

# remove data with NAS in relevant variables
idx_noNA <- which(rowSums(is.na(all_final[, !names(all_final) %in% c("imgincl_t1w_include","mrif_score","rel_family_id")])) == 0)

all_final_noNA <- all_final[idx_noNA, ]

# change the data type
str(all_final_noNA)

all_final_noNA <- all_final_noNA %>% mutate_at(c("sex","site","race_ethnicity","parental_education"), as.factor)
str(all_final_noNA)
colSums(is.na(all_final_noNA))

# adjust the level of parental education
levels(all_final_noNA$parental_education) <- list("low" = 1:12,"medium" = 13:17,"high" = 18:21)
str(all_final_noNA)
dim(all_final_noNA)

saveRDS(all_final_noNA,"all_final_noNA_item.rds")

names(all_final_noNA)
########################################
#train test split: site
#########################################
# first split the site 10 times
sites <- unique(all_final_noNA$site)
# there are 22 sites, so 22*0.8=17.6 for the training sites
train_sites <- lapply(1:10, function(i) sample(sites, size = 17, replace = FALSE))

# then split of the ids
all_final_train_item <- lapply(1:10, function(i) all_final_noNA %>% filter(site %in% train_sites[[i]]))
all_final_test_item <- lapply(1:10, function(i) all_final_noNA %>% filter(!site %in% train_sites[[i]]))

lapply(all_final_train_item, nrow)
class(all_final_train[[1]])
saveRDS(all_final_train_item,"all_final_train_item.rds")
saveRDS(all_final_test_item,"all_final_test_item.rds")















