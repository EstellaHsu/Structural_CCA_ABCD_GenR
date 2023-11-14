################################################################################################################################
################################ further validation of the brain canonical variates in GenR ####################################
################################################################################################################################


##############################################
############## school achievement @13 ########
##############################################

###### read the data
# extract school achievement
cito <- read.csv("GR1093-F_School_13052020.csv",header = TRUE)
names(cito)[1] <- "IDC"
cito1 <- cito[,c("IDC","F0200293_cleaned","F0200293_cleaned_CITOscore")]

# extract CV scores
CV_scores <- scale(genr_brain_str) %*% brain_mean # brain_mean is the averaged canonical loadings from 10 train-test splits in ABCD
CV_scores <- as.data.frame(CV_scores)
names(CV_scores) <- c(paste0("CV",1:3),"IDC")

# merge all the data
df_cv_school <- merge(cito1, CV_scores, by="IDC")
df_cv_school <- merge(df_cv_school, demo[,c("IDC","age","gender","ethnicity","maternal_edu")])
dim(df_cv_school)
names(df_cv_school)[3] <- "school_achievement"

cito_cv1 <- lm(scale(school_achievement) ~  scale(CV1) + age + gender + ethnicity + maternal_edu, df_cv_school)
summary(cito_cv1)
confint(cito_cv1)


##############################################
################# cognition @13 ##############
##############################################
wisc13 <- read.csv("CHILDWISC13_16082021.csv", header = TRUE)
names(wisc13)[1] <- "IDC"

wisc_all <- merge(df_cv_school, wisc13[,c("IDC","WISC13_Voc_Tscore","WISC13_MR_Tscore",
                                   "WISC13_DS_Tscore","WISC13_CD_Tscore","WISC13_FSIQ")], by="IDC")

colSums(is.na(wisc_all))
idx_na_wisc <- which(rowSums(is.na(wisc_all)) != 0)
wisc_all_0 <- wisc_all[-idx_na_wisc, ]
colSums(is.na(wisc_all_0))
names(wisc_all_0)


voc_cv1 <- lm(scale(WISC13_Voc_Tscore) ~ scale(CV1) + age + gender + ethnicity + maternal_edu, wisc_all)
summary(voc_cv1)
confint(voc_cv1)

mr_cv1 <- lm(scale(WISC13_MR_Tscore) ~ scale(CV1) + age + gender + ethnicity + maternal_edu, wisc_all)
summary(mr_cv1)
confint(mr_cv1)

ds_cv1 <- lm(scale(WISC13_DS_Tscore) ~ scale(CV1) + age + gender + ethnicity + maternal_edu, wisc_all)
summary(ds_cv1)
confint(ds_grey)

cd_cv1 <- lm(scale(WISC13_CD_Tscore) ~ scale(CV1) + age + gender + ethnicity + maternal_edu, wisc_all)
summary(cd_cv1)
confint(cd_cv1)

fsIQ_cv1 <- lm(scale(WISC13_FSIQ) ~ scale(CV1) + age + gender + ethnicity + maternal_edu, wisc_all)
summary(fsIQ_cv1)
confint(fsIQ_cv1)

########################################################
################# ADHD medication use @10 ##############
########################################################

adhd <- read.csv("ADHDmed.csv")
names(adhd)[1] <- "IDC"
df_adhd <- merge(df_cv_school, adhd, by="IDC")
dim(df_adhd)
str(df_adhd)

df_adhd[df_adhd == 999] <- NA
df_adhd$ADHDmed <- as.factor(df_adhd$ADHDmed)

summary(glm(ADHDmed ~ scale(CV1) + age + gender + ethnicity + maternal_edu,family=binomial(), df_adhd))
confint(glm(ADHDmed ~ scale(CV1) + age + gender + ethnicity + maternal_edu,family=binomial(), df_adhd))


########################################################
################# ADHD polygenic risk scores ###########
########################################################
# PGS for ADHD
prs_adhd <- fread("ADHD_GENR.prs.all_score")
names(prs_adhd)[2] <- "IDC"

# PCs
pcs <- fread("PCA_SelectionGWAv3_Cauc_May2014.csv")

df1 <- merge(df_cv_school,prs_adhd[,c("IDC","Pt_1")], by="IDC")
df_prs_adhd <- merge(df1,pcs[,c("IDC","C1","C2","C3","C4","C5")], by="IDC")
names(df2)[9] <- "prs_adhd"

prs_adhd <- lm(scale(CV1) ~ scale(prs_adhd) + age + gender + C1 + C2 + C3 + C4 + C5, df_prs_adhd)
summary(prs_adhd)
confint(prs_adhd)
