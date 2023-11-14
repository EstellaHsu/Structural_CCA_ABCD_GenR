##############################################
########## school achievement prediction #####
##############################################

###### read the data
setwd("V:/medewerkers/****** Xu, B/PhD projects/Structural_CCA/")
# extract school achievement
cito <- read.csv("GR1093-F_School_13052020.csv",header = TRUE)
names(cito)[1] <- "IDC"

cito1 <- cito[,c("IDC","F0200293_cleaned","F0200293_cleaned_CITOscore")]

#extract the ids
demo <- readRDS("all_final_genr_tbv.rds")
id <- demo$IDC
dim(demo)
names(demo)

# extract CV scores
CV_scores <- readRDS("CV_genr_total_31.01.2023.rds")
CV_scores <- as.data.frame(CV_scores)
dim(CV_scores)
CV_scores$IDC <- id
names(CV_scores) <- c(paste0("CV",1:3),"IDC")

# merge all the data
df_cv_school <- merge(cito1, CV_scores, by="IDC")
df_cv_school <- merge(df_cv_school,global1, by="IDC")
names(df_cv_school)
df_cv_school <- merge(df_cv_school, demo[,c("IDC","age","gender","ethnicity","maternal_edu","genr_tbv_f09")])
dim(df_cv_school)
names(df_cv_school)[3] <- "school_achievement"

cito_cv1 <- lm(scale(school_achievement) ~  scale(CV1) + age + gender + ethnicity + maternal_edu, df_cv_school)
summary(cito_cv1)
cito_eTIV <- lm(scale(school_achievement) ~  scale(eTIV_f09) + age + gender + ethnicity + maternal_edu, df_cv_school)
summary(cito_eTIV)
cito_totalgrey <- lm(scale(school_achievement) ~  scale(genr_tbv_f09) + age + gender + ethnicity + maternal_edu, df_cv_school)
summary(cito_totalgrey)

data <- data.frame(pred=predict(cito_cv3),
                   actual=df_cv_school[!is.na(df_cv_school$school_achievement),"school_achievement"])
mean((data$pred - data$actual)^2)


confint(cito_cv1)






