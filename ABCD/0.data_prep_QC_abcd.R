###########################################################################################################################
######################## ABCD: data preparation based on brain data QC, covariates, and CBCL ##############################
###########################################################################################################################

# first read all the QC and covariates data
# a small function to rename all the ids to match the ids of the brain data in our server
rename_idc <- function(x){ sapply(as.character(x$src_subject_id), 
                           function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))}

#########################################
#demographic info
#########################################
participants <- read.delim(file = 'abcd_lpds01.txt')
names(participants)
participants$idc <- rename_idc(participants)

# only 1 year follow up had the parental education data
prnt_ed <- participants[-1,c("idc","demo_prnt_ed_v2_l","eventname")]
prnt_ed1 <- prnt_ed %>% filter(as.character(eventname) == "1_year_follow_up_y_arm_1")

# other data we need baseline (age 10)
participants1 <- participants %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
# extract ids, sex, and age data
participants2 <- participants1[, c("idc","sex","interview_age")]

# merge the useful variables
demo_info <- merge(prnt_ed1[,c("idc","demo_prnt_ed_v2_l")],participants2,by="idc")

#########################################
#Site information
#########################################
site <- read.delim(file = 'abcd_lt01.txt')
names(site)                                
site$idc <- rename_idc(site)
site1 <- site %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
# extract the study site information
site_info <- site1[, c("idc","site_id_l")]


#########################################
#twins & siblings information
#########################################
twin <- read.delim("acspsw03.txt")
names(twin)

twin$idc <- rename_idc(twin)                  
twin1 <- twin %>% filter(eventname == "baseline_year_1_arm_1")
# extract the race & ethnicity, and family ids (used for exlude twins/siblings/triples)           
twin_info <- twin1[,c("idc","race_ethnicity","rel_family_id")]


#########################################
#brain data QC
#########################################
                   
qc_str <- read.delim("abcd_imgincl01.txt")
names(qc_str)

qc_str$idc <- rename_idc(qc_str)
qc1 <- qc_str %>% filter(eventname == "baseline_year_1_arm_1")
# extract the recommended QC variable
strqc_info1 <- qc1[,c("idc","imgincl_t1w_include")]


#########################################
#brain data incidental findings
#########################################
setwd("/Users/estella/Desktop/ABCD_download")
qc_findings <- read.delim("abcd_mrfindings02.txt")
names(qc_findings)

qc_findings$idc <- rename_idc(qc_findings)
qc2 <- qc_findings %>% filter(eventname == "baseline_year_1_arm_1")

#extract the incidental finding scroes
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

# 1. ids who have the brain structral data on our server
id <- readRDS("id_has_strdata.rds")
all_info <- all_info %>% filter(idc %in% id)
dim(all_info)

# 1. inclusion criteria: based on the inclusion critaria of ABCD and remove the children with clinical relevant incidental findings
all_info1 <- all_info %>% filter(imgincl_t1w_include == 1 & mrif_score < 3)
dim(all_info1)

# 2. remove twin/siblings/triples
all_info2 <- all_info1[!duplicated(all_info1$rel_family_id), ]
dim(all_info2)


                                  
#########################################
#merge CBCL: behavioral assessment
#########################################
cbcl <- read.delim("abcd_cbcls01.txt")
cbcl$idc <- rename_idc(cbcl)

cbcl1 <- cbcl %>% filter(as.character(eventname) == "baseline_year_1_arm_1")
names(cbcl1)
                                  
                                  
######## extract syndrome scale scores, internalizing, externalizing and total problems raw scoress
cbcl_syn <- cbcl1[, c("idc","cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                      "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                      "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                      "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r")]
dim(cbcl_syn)


###### merge all the data: cbcl & post-qc data

all_cbcl <- merge(all_info2, cbcl_syn,by="idc")
str(all_cbcl)

# the structrue of the original data is all characters, so we transformed some of them to numeric data
all_cbcl1 <- all_cbcl %>% mutate_at(c("interview_age","cbcl_scr_syn_anxdep_r","cbcl_scr_syn_withdep_r","cbcl_scr_syn_somatic_r",
                         "cbcl_scr_syn_social_r","cbcl_scr_syn_thought_r","cbcl_scr_syn_attention_r",
                         "cbcl_scr_syn_rulebreak_r","cbcl_scr_syn_aggressive_r","cbcl_scr_syn_internal_r",
                         "cbcl_scr_syn_external_r","cbcl_scr_syn_totprob_r"),as.numeric) 

colSums(is.na(all_cbcl1))

# there are many NAs in the final data in parental education, sex, interview age: over 400
# so I tried to read the ABCC participants data to find more demographic information (baseline at age 10)
dd_rsmri <- read.table(file = 'participants.tsv', sep = '\t', header = TRUE)
names(dd_rsmri)[1] <- "idc"
# use the parental edu, sex, and age from this dataset. 
rsmri <- dd_rsmri[, c("idc","sex","age","parental_education","site")]

                                  
                                  
# fill the missingness
na_demo <- all_cbcl1[is.na(all_cbcl1$demo_prnt_ed_v2_l), ]
dim(na_demo)
na_fill <- merge(rsmri, na_demo[,!names(na_demo) %in% c("sex","interview_age","demo_prnt_ed_v2_l")], by=c("idc","site"))
dim(na_fill)

# replace the missingness 
all_final0 <- all_cbcl1 %>% 
  mutate(demo_prnt_ed_v2_l=replace(demo_prnt_ed_v2_l, is.na(demo_prnt_ed_v2_l), na_fill$parental_education)) %>%
  mutate(interview_age=replace(interview_age, is.na(interview_age), na_fill$age)) %>%
  mutate(sex=replace(sex, is.na(sex), na_fill$sex))

dim(all_final0)
all_final0$sex[all_final0$sex == 1] <- "M"
all_final0$sex[all_final0$sex == 2] <- "F"

names(all_final0)[6:7] <- c("site","parental_education")

all_final0[all_final0 == 777] <- NA
all_final0[all_final0 == 888] <- NA
all_final0[all_final0 == 999] <- NA
all_final0[all_final0 == ""] <- NA

colSums(is.na(all_final0))

# remove data with NAS in relevant variables
idx_noNA <- which(rowSums(is.na(all_final0[, !names(all_final0) %in% c("imgincl_t1w_include","mrif_score","rel_family_id")])) == 0)

all_final_noNA <- all_final0[idx_noNA, ]
dim(all_final_noNA)

# change the data type
str(all_final_noNA)

all_final_noNA <- all_final_noNA %>% mutate_at(c("sex","site","race_ethnicity","parental_education"), as.factor)
str(all_final_noNA)

# adjust the level of parental education
levels(all_final_noNA$parental_education) <- list("low" = 1:12,"medium" = 13:17,"high" = 18:21)
str(all_final_noNA)
dim(all_final_noNA)

saveRDS(all_final_noNA,"all_final_noNA_abcd.rds")


########################################
#train test split based on study sites
#########################################

set.seed(519)
# first split the site 10 times
sites <- unique(all_final_noNA$site)
# there are 21 sites, so 21*0.8=16.8 for the training sites
train_sites <- lapply(1:10, function(i) sample(sites, size = 17, replace = FALSE))

# then split of the ids
all_final_train <- lapply(1:10, function(i) all_final_noNA %>% filter(site %in% train_sites[[i]]))
all_final_test <- lapply(1:10, function(i) all_final_noNA %>% filter(!site %in% train_sites[[i]]))


saveRDS(all_final_train,"all_final_train.rds")
saveRDS(all_final_test,"all_final_test.rds")


                         
                        
