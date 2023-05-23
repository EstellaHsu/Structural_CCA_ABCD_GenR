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
all_abcd <- readRDS("train_test_split.rds")

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
  
  brain_train <- str_train[[i]]$brain_train
  cbcl_train <- str_train[[i]]$cbcl_train
  brain_test <- str_test[[i]]$brain_test
  cbcl_test <- str_test[[i]]$cbcl_test
  brain_pen <- grid.abcd.str.all[[i]][[1]][1]
  cbcl_pen <- grid.abcd.str.all[[i]][[1]][2]
  
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  perm_abcd_train <- permutation_test(cbcl_train, brain_train, nperm=999, cbcl_pen, brain_pen, 8, res.abcd$cors)
  
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd, 8)
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=1999,res.abcd,abcd.test)

  cor.abcdTogenr <- test_project_weights(brain_genr_str, cbcl_genr_str, res.abcd,8)
  perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr_str,brain_genr_str,nperm=1999, 
                                                    res.abcd,abs(cor.abcdTogenr))
  
  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr_total$pval.perm))
})


saveRDS(str_train_test_abcd,"str_train_test_abcd_31.01.2023.rds")


#################################################
############## traditional CCA ##################
#################################################

## I just did the traditional CCA in the first train-test split (it can be done in any split)


brain_train <- str_train[[1]]$brain_train
cbcl_train <- str_train[[1]]$cbcl_train
brain_test <- str_test[[1]]$brain_test
cbcl_test <- str_test[[1]]$cbcl_test

## the traditional CCA model
abcd_train_cca <- candisc::cancor(brain_train, cbcl_train)


## calculate the canonical correlations in the test set of ABCD
std_brain <- scale(brain_test) %*% abcd_train_cca$coef$X
std_cbcl <- scale(cbcl_test) %*% abcd_train_cca$coef$Y

cor.res.test.abcd <- diag(cor(std_brain, std_cbcl))

## test the significance using permutation tests
# number of permutations
n.perm = 1999

# shffule the row the behavioral data and build a new list of cbcl data 
shuffle_idx <- sapply(1:n.perm, function (x){permute::shuffle(1:nrow(cbcl_test))})
cbcl_perm <- lapply(1:n.perm, function(i) {cbcl_test[shuffle_idx[, i], ]})

# build the null distribution of the canonical correlations
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
perm_cor <- foreach::foreach(i = seq_along(cbcl_perm)) %dopar% {
  std_brain <- scale(brain_test) %*% abcd_train_cca$coef$X
  std_cbcl <- scale(cbcl_perm[[i]]) %*% abcd_train_cca$coef$Y
  cors <- diag(cor(std_brain, std_cbcl))
}
stopCluster(cl)

cor.perm <- do.call(rbind, perm_cor)
# calculate the p values 
pval.perm <- sapply(1:6, function(x){length(which(abs(cor.perm[, x]) >= abs(cor.res.test.abcd[x])))/(n.perm+1)})


## test the correlations in generation R
std_brain <- scale(brain_genr_str) %*% abcd_train_cca$coef$X
std_cbcl <- scale(cbcl_genr_str) %*% abcd_train_cca$coef$Y
cor.res.genr <- diag(cor(std_brain, std_cbcl))




########################################################
############## reorder and average across splits  ######
########################################################

str_train_test_abcd <- readRDS("str_train_test_abcd.rds")
brain_mean <- readRDS("brain_mean.rds")
cbcl_mean <- readRDS("cbcl_mean.rds")

# the reference split
load_std <- res.abcd$v # reorder based on brain loadings
# we only care about the first 4 CVs based on the covariance explained plot. 
# As the CBCL loadings are high similar, if we have more CVs, the brain might be chaotic

new_cv_reorder  <- lapply(seq_along(str_train_test_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:5) {
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:5, function(z) {cor(load_std[, i], str_train_test_abcd[[1]]$abcd.train$v[, z])})
    if(max(abs(cor_cv)) < 0.5){ # sometimes the correlation is too low
      idx  <- i
    } else {
      idx <- which(abs(cor_cv) == max(abs(cor_cv)))
    }
    
    if(!idx %in% idx_v) {
      idx_v <- c(idx_v, idx) # sometimes it will overlap
    } else {
      idx_v <- c(idx_v, i) 
    }
  }
  idx_v[duplicated(idx_v)] <- c(1:5)[!c(1:5) %in% idx_v]
  
  loadMat <- str_train_test_abcd[[j]]$abcd.train$v[, idx_v]
  return(list(loadMat=loadMat, idx_reorder=idx_v))
  
})


lapply(new_cv_reorder, function(x) {x$idx_reorder})
c <- lapply(new_cv_reorder, function(x) {x$loadMat})
cbcl_mean <- do.call(cbind, lapply(1:5, 
                                function(i) {
                                  rowMeans(abs(do.call(cbind, 
                                                       lapply(c, function(x) {x[, i]})
                                  )))}))

rownames(cbcl_mean)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")

corrplot(t(abs(cbcl_mean))[1:3,], method="color", 
         addCoef.col = "black", tl.srt =45, 
         tl.col = "black", tl.cex = 2, number.cex=1)
##### brain average
b <- lapply(seq_along(new_cv_reorder), function(x) {
  loadMat <- str_train_test_abcd[[x]]$abcd.train$u[, new_cv_reorder[[x]]$idx_reorder]})

brain_mean <- do.call(cbind, lapply(1:5, 
                                    function(i) {
                                      rowMeans(abs(do.call(cbind, lapply(b, function(x) {x[, i]})
                                      )))}))

######### as the signn of the loadings could be changed across splits
######### we aim to match the signs with the reference plit
signC <- sign(res.abcd$v[,1:3]) #cbcl
signB <- sign(res.abcd$u[,1:3]) # brain

# as it is very complicated to match every sign of the loadings, we only match the most important one
signB.cv1 <- signB[which.max(abs(res.abcd$u[,1])),1]
signB.cv2 <- signB[which.max(abs(res.abcd$u[,2])),2]
signB.cv3 <- signB[which.max(abs(res.abcd$u[,3])),3]

# the signs are relative within one canonical variate, if we set the loadings for CBCL as positive
brain_mean_sign <- brain_mean * -1


rownames(brain_mean) <- colnames(str_train[[8]]$brain_train)
brain_check <- round(brain_mean, 3)
quantile(brain_mean)
rownames(brain_check) <- colnames(str_train[[8]]$brain_train)
brain_check[-c(69:102),1:3]
brain_check[brain_check < 0.08351686] <- 0
saveRDS(cbcl_mean, "cbcl_mean.rds")
saveRDS(brain_mean, "brain_mean.rds")
brain_mean2 <- readRDS("brain_mean.rds")





cv_abcd <- scale(str_test[[8]]$brain_test) %*% brain_mean
cbcl_cv_abcd <- scale(str_test[[8]]$cbcl_test) %*% cbcl_mean
cv_abcd <- round(cv_abcd,3)


str_test[[8]]$brain_test
cor(cbcl_cv_abcd[,1], cbcl_cv_abcd[,3])
hist(cv_abcd[,1])
cor(brain_check[,])


plot(cv_abcd[,3], cbcl_cv_abcd[,3])

cor(cbcl_cv_abcd[,2],scale(str_test[[8]]$cbcl_test))

########################################################
############## train-test correlations ################
########################################################


# reorder the correlations
traincor <- lapply(str_train_test_abcd, function(x) {x$abcd.train$cors})
traincor <- do.call(rbind,traincor)
traincor_reorder <- lapply(1:10, function(i) {traincor[i, new_cv_reorder[[i]]$idx_reorder]})
traincor_reorder <- do.call(rbind,traincor_reorder)
rownames(traincor_reorder) <- paste0("split",1:10)
colnames(traincor_reorder) <- paste0("CV",1:5)

meantrain <- colMeans(traincor_reorder)
sdtrain <- apply(traincor_reorder, 2, sd)


testcor <- lapply(str_train_test_abcd, function(x) {x$abcd.test})
testcor <- do.call(rbind,testcor)
testcor_reorder <- lapply(1:10, function(i) {testcor[i, new_cv_reorder[[i]]$idx_reorder]})
testcor_reorder <- do.call(rbind,testcor_reorder)
rownames(testcor_reorder) <- paste0("split",1:10)
colnames(testcor_reorder) <- paste0("CV",1:5)
meantest <- colMeans(testcor_reorder)
sdtest <- apply(testcor_reorder, 2, sd)

fig_10splits_syn <- data.frame(Train_Test = rep(c("Training set","Test set"),each=3), meancor = c(meantrain[1:3], meantest[1:3]),
                               sdcor = c(sdtrain[1:3], sdtest[1:3]), CV=rep(paste0("CV",1:3),2))


p <- ggplot(fig_10splits_syn, aes(x=CV, y=meancor, color=Train_Test, group=Train_Test)) + 
  scale_color_manual(name="Training Test Set", values = c("#7D3C98", "#52BE80")) + 
  geom_errorbar(aes(ymin=meancor-sdcor, ymax=meancor+sdcor), width=.2, size=3) +
  geom_point(size=8) + ylim(c(0, 0.2)) + 
  theme_bw() + labs(y = "canonical correlations", x = "canonical variates", size=20) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.title = element_text(size=20,face="bold"),
        legend.text = element_text(size=20),
        legend.position = c(0.8, 0.86))
p

ggsave("train_test_variance_03.03.2023.pdf")

tab <- round(brain_mean[, 1:3],3)
rownames(tab)
quantile(tab[103:109, ]) # 0.089, 0.088, 0.086

########################################################
############## CBCL radar Plot #########################
########################################################

temp <- cbcl_mean
df_temp1 <- as.data.frame(temp)
df_temp1 <- t(df_temp1)
df_temp2 <- rbind(rep(1,8), df_temp1)
colnames(df_temp2) <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule breaking","Aggression")

df_temp3 <- as.data.frame(cbind(group=c("1","CV1", "CV2", "CV3", "CV4", "CV5"),df_temp2))
df_temp3[,2:9] <- lapply(df_temp3[,2:9], function(x){as.numeric(as.character(x))})
df_temp3 <- df_temp3[1:4,]

#pdf("radarfigure_train_abcd.pdf", width=12, height=8)
colors_border_6 <- c("#fcfcfb","#FFDF6F","#B6DEE7","#D65353")
colors_in_6 <- c("#fcfcfb","#FFDF6F","#B6DEE7","#D65353")

p2 <- ggradar2(df_temp3, base.size=5, webtype = 'lux', grid.min = 0,
               grid.max = 1,label.gridline.mid = TRUE, group.colours = colors_border_6,
               group.fill.colours = colors_in_6,label.centre.y=FALSE,
               gridline.mid.colour="grey", grid.label.size = 0,
               gridline.max.linetype = "solid",polygonfill.transparency=0.5,
               background.circle.transparency=0.1,axis.label.size=5.5,
               group.line.width=1.5,group.point.size=3,plot.legend = TRUE) 
p2

ggsave("radarfigure_str_abcd.pdf", width = 12, height = 6)
dev.off()


########################################################
############## Covariance explained  ###################
########################################################
vardf <- VarianceExplain(brain_test, cbcl_test, res.abcd, 8) 
colnames(vardf) <- c("Canonical_Variates", "Covariance_Explained")
vardf$`Canonical_Variates` <- paste0("CV",1:8)
vardf$Canonical_Variates <- as.factor(vardf$Canonical_Variates)
ggplot(vardf, aes(x=Canonical_Variates, y=Covariance_Explained, group=Canonical_Variates))+
  geom_point(aes(shape=Canonical_Variates,color=Canonical_Variates), size = 7, alpha = 0.9) + 
  theme_bw() + theme(legend.position="none",
                     axis.text.x = element_text(size = 12)) + 
  scale_shape_manual(values=c(17,16,15,8,8,8,8,8)) +
  scale_color_manual(values=c("#FFDF6F","#B6DEE7","#D65353","grey","grey","grey","grey","grey")) 

ggsave("covariance_explained_testset.pdf",width=8,height=4)


### permutation test: an example
perm_abcd7 <- permutation_test(cbcl_train, brain_train, nperm=999, 0.5, 0.5, 8, res.abcd.7$cors)


########################################################
############## permutation test visualization  #########
########################################################

# direct mapping
cor.abcd.test <- test_project_weights(brain_test, cbcl_test, res.abcd, 8)
cor.abcd.test
perm_abcdtest <- permutation_test_testset(cbcl_test, brain_test,nperm=999, 
                                           res.abcd,cor.abcd.test)
#"#A0D568","#F7EA48","#4FC1E8"
# direct mapping
cor.abcdTogenr <- test_project_weights(brain_genr_str, cbcl_genr_str, res.abcd, 6)
cor.abcdTogenr

perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr_str,brain_genr_str,nperm=999, 
                                                  res.abcd,abs(cor.abcdTogenr))

p.adjust(perm_abcdtest$pval.perm)


df.cor.2 <- data.frame(cor = perm_abcdtest$cor.perm[, 2])
df.cor.1 <- data.frame(cor = perm_abcdtest$cor.perm[, 1])
df.cor.3 <- data.frame(cor = perm_abcdtest$cor.perm[, 3])

ggplot(df.cor.3, aes(x=cor)) +  # Apply nrow function
  geom_histogram(binwidth=0.008, fill = "#D65353", color = "black", alpha = 0.8)+
  geom_vline(xintercept=0.113, linetype="dashed", color = "red", size = 1) +
  ggtitle("Permutation test of the 1st canonical correlation") + 
  labs(x = "canonical correlations") + 
  #annotate(geom="text", x=0.0, y=200, label="r = 0.13", size = , color="black") + 
  annotate(geom="text", x=0.06, y=120, label="P(FDR) < 0.001", 
           fontface = 'italic', size = 6, color="black") + theme_bw()+
  theme(
    #panel.background = element_rect(fill = "#F5F7FA",
                                    #colour = "#F5F7FA",
                                    #size = 0.5, linetype = "solid"),
    #panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                    #colour = "white"), 
    #panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                    #colour = "white"),
    axis.text=element_text(size=20),
    axis.title=element_text(size=20,face="bold"),
    legend.title = element_text(size=20,face="bold"),
    legend.text = element_text(size=20),
    plot.title = element_text(size=22, face="bold")
  )
ggsave("perm_abcd_3.pdf",width=10,height=5)





########################################################
############## brain figures ###########################
########################################################
brain_fig <- brain_mean[,1:3]
quantile(brain_fig[35:68, ])

# choose the 75% quantile
brain_fig <- round(brain_fig, 3)

rownames(brain_fig) <- colnames(str_train[[8]]$brain_train)
colnames(brain_fig) <- paste0("CV",1:3)
#brain_fig <- brain_fig[-which(rowSums(brain_fig) == 0),]

library(ggseg)
library(ggseg3d)
ggseg(mapping=aes(fill = region), colour="black")

rownames(brain_fig)
# create a figure data frame of surface areas
fig_surf <- brain_fig[35:68, ]
quantile(fig_surf)
fig_surf[fig_surf < 0.11325] <- 0
surf <- gsub("_surfarea","",rownames(fig_surf))
# some of the names are not correct
surf <- c("bankssts","caudal anterior cingulate", "caudal middle frontal","cuneus","entorhinal","fusiform",
          "inferior parietal","inferior temporal","isthmus cingulate","lateral occipital","lateral orbitofrontal","lingual",
          "medial orbitofrontal","middle temporal", "parahippocampal","paracentral","pars opercularis",
          "pars orbitalis","pars triangularis","pericalcarine","postcentral","posterior cingulate","precentral",
          "precuneus","rostral anterior cingulate","rostral middle frontal","superior frontal","superior parietal",
          "superior temporal","supramarginal","frontal pole","temporal pole","transverse temporal","insula")
length(surf)
fig_surf <- as.data.frame(fig_surf)
fig_surf$region <- surf
names(fig_surf)

######## 2D figures
ggseg(.data=fig_surf[,c("region","CV3")], mapping=aes(fill=CV3),colour="grey") +
  labs(title="CV3 Surface areas", fill="loadings") +
  scale_fill_gradient(low="white",high="#D65353",na.value="#F2F3F4")

ggsave("CV3_surfarea.pdf",width=10,height=4)


######## 3D figures
ggseg3d(.data = fig_surf[,c(3,4)], 
        atlas = dk_3d, 
        na.alpha= .5, 
        colour = "CV3",
        palette = c("#F8F9F9" = 0, "#D65353" = .2, "#884EA0" = 0.3)) %>% 
  remove_axes()


# create a figure data frame of cortical volumes
fig_vol <- brain_fig[1:34, ]
quantile(fig_vol)
fig_vol[fig_vol < 0.10875] <- 0
fig_vol <- as.data.frame(fig_vol)
fig_vol$region <- surf
names(fig_vol)

# CV1
ggseg(.data=fig_vol[,c("region","CV3")], mapping=aes(fill=CV3),colour="grey") +
  labs(title="CV3 Cortical volumes", fill="loadings") +
  scale_fill_gradient(low="white",high="#D65353",na.value="#F2F3F4")

ggsave("CV3_vol.pdf",width=10,height=4)

#low="#AED6F1",high="#2980B9"
#low="#F7DC6F",high="#D35400"

ggseg(atlas = "aseg",mapping=aes(fill = region), colour="black")

# create a figure data frame of surface areas
fig_sub <- brain_fig[47:51, ]
sub <- gsub("_vol","",rownames(fig_sub))
sub <- c("hippocampus","amygdala","accumbens","thalamus","putamen")
rownames(fig_surf) <- surf


########################################################
############## brain figures: subcortical ##############
########################################################

plot(aseg)
rownames(brain_fig)
# create a figure data frame of cortical volumes
fig_sub <- brain_fig[103:109, ]
fig_sub[fig_sub < 0.08] <- 0
quantile(fig_sub)

# CV1
cv_loadings <- fig_sub[,3]
subcor <- gsub("_vol","",rownames(fig_sub))
subcor <- c("hippocampus", "amygdala","accumbens area","caudate","thalamus proper","putamen","pallidum")
fig_str <- data.frame(region = subcor, loadings=fig_sub[,3])
#fig_str <- fig_str[fig_str$loadings != 0,]

ggseg(.data=fig_str, atlas=aseg,mapping=aes(fill=loadings),colour="white",view = "coronal", label=TRUE) +
  labs(title="CV1 Subcortical volumes", fill="loadings") +
  scale_fill_gradient(low="#FFDF6F",high="#F39C12")

ggsave("CV1_sub.pdf",width=4,height=4)


######## 3d subcortical

# CV1
cv_loadings <- fig_sub[,3]

fig_aseg <-  aseg_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  dplyr::select(label) %>% 
  filter(grepl("Hippocampus|Amygdala|Accumbens-area|Caudate|Thalamus-Proper|Putamen|Pallidum", 
                label)) %>% 
  mutate(loadings = rep(cv_loadings, 2)) 

p <- ggseg3d(.data = fig_aseg, atlas = aseg_3d, 
        colour = "loadings", text = "loadings", 
        na.alpha= .5,
        palette = c("white" = 0, "Salmon" = .075, "#CB4335" = 0.15)) %>% 
  remove_axes() %>%
  add_glassbrain(colour = "white",
                 opacity = 0.2)
p

########################################################
############## ABCD wisc added: what's the difference ##
########################################################
setwd("~/Desktop/Structural_cca/")

allinfo <- readRDS("all_final_noNA_update.rds")
id <- allinfo$idc

############# read the pearson data
setwd("/Users/estella/Desktop/ABCD_download/abcd_pearson")
pearson <- read.delim("abcd_ps01.txt")
names(pearson)

pearson$idc <- sapply(as.character(pearson$src_subject_id), 
                                 function(x) paste0("sub-NDAR", strsplit(x, "_")[[1]][2]))

pearson <- pearson %>% filter(eventname == "baseline_year_1_arm_1")

pearson_abcd <- pearson[pearson$idc %in% id, c("idc","pea_wiscv_tss")] %>% mutate(wisc=as.numeric(pea_wiscv_tss))
dim(pearson_abcd)
pearson_abcd[pearson_abcd == 888] <- NA
pearson_abcd[pearson_abcd == 999] <- NA
colSums(is.na(pearson_abcd)) # no missing

############# do the grid search
setwd("~/Desktop/Structural_cca/")
id_train <- readRDS("all_final_train.rds")
id_train <- lapply(id_train, function(x) {x$idc})
id_test <- readRDS("all_final_test.rds")
id_test <- lapply(id_test, function(x) {x$idc})
# check whether the order of ids is identical
#identical(as.character(id_train[[1]]),wisc$idc)
# not identical, so we need to reorder
grid.abcd.str.wisc <- lapply(1:10, function(i) {
  
  structure_resample <- str_train[[i]]
  wisc <- pearson_abcd %>% filter(idc %in% id_train[[i]]) %>% mutate(wisc=as.numeric(pea_wiscv_tss))
  wisc <- wisc[order(wisc$idc),]
  
  cbcl_train <- structure_resample$cbcl_train
  cbcl_train$wisc <- wisc$wisc
  cbcl_wisc_train <- cbcl_train
  print(nrow(cbcl_wisc_train))
  id_na <- which(is.na(cbcl_wisc_train$wisc))
  cbcl_wisc_train <- cbcl_wisc_train[-id_na,]
  print(nrow(cbcl_wisc_train))
  brain_train <- structure_resample$brain_train
  brain_train <- brain_train[-id_na,]
  print(nrow(brain_train))
  
  cv.resample <- CV_sampling(cbcl_wisc_train, brain_train, 100, 0.8)
  
  x_pen <- seq(0.1, 1.0, length.out = 10) # brain
  y_pen <- seq(0.1, 1.0, length.out = 10)
  
  grid.abcd <- grid.search.cor.Testset(cv.resample$brain_train, cv.resample$cbcl_train,
                                       cv.resample$brain_test, cv.resample$cbcl_test,
                                       x_pen, y_pen, nsample=100)
})

saveRDS(grid.abcd.str.wisc, "grid.abcd.str.wisc.rds")

############### fit sCCA model
#cbcl_pen <- c(0.6,0.7,0.6,0.6,0.7,0.7,0.7,0.6,0.7,0.7)
#brain_pen <- c(0.7,0.5,0.7,0.7,0.4,0.6,0.5,0.6,0.6,0.3)

str_train_test_wisc_abcd <- lapply(1:10, function(i) {
  
  structure_resample <- str_train[[i]]
  wisc <- pearson_abcd %>% filter(idc %in% id_train[[i]]) %>% mutate(wisc=as.numeric(pea_wiscv_tss))
  wisc <- wisc[order(wisc$idc),]
  
  cbcl_train <- structure_resample$cbcl_train
  cbcl_train$wisc <- wisc$wisc
  cbcl_wisc_train <- cbcl_train
  print(nrow(cbcl_wisc_train))
  id_na <- which(is.na(cbcl_wisc_train$wisc))
  cbcl_wisc_train <- cbcl_wisc_train[-id_na,]
  print(nrow(cbcl_wisc_train))
  brain_train <- structure_resample$brain_train
  brain_train <- brain_train[-id_na,]
  print(nrow(brain_train))
  
  
  structure_resample_test <- str_test[[i]]
  wisc0 <- pearson_abcd %>% filter(idc %in% id_test[[i]]) %>% mutate(wisc=as.numeric(pea_wiscv_tss))
  wisc0 <- wisc0[order(wisc0$idc),]
  
  cbcl_test <- structure_resample_test$cbcl_test
  cbcl_test$wisc <- wisc0$wisc
  cbcl_wisc_test <- cbcl_test
  print(nrow(cbcl_wisc_test))
  id_na <- which(is.na(cbcl_wisc_test$wisc))
  cbcl_wisc_test <- cbcl_wisc_test[-id_na,]
  print(nrow(cbcl_wisc_test))
  brain_test <- structure_resample_test$brain_test
  brain_test <- brain_test[-id_na,]
  print(nrow(brain_test))
  
  brain_pen <- grid.abcd.str.wisc[[i]][[1]][[1]]
  cbcl_pen <- grid.abcd.str.wisc[[i]][[1]][[2]]
  
  res.abcd <- CCA(x=brain_train, z=cbcl_wisc_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  perm_abcd_train <- permutation_test(cbcl_wisc_train, brain_train, nperm=999, cbcl_pen, brain_pen, 8, res.abcd$cors)
  
  abcd.test <- test_project_weights(brain_test,cbcl_wisc_test,res.abcd, 8)
  perm_abcd_test <- permutation_test_testset(cbcl_wisc_test,brain_test,nperm=999,res.abcd,abcd.test)
  
  cor.abcdTogenr <- test_project_weights(brain_genr_str, cbcl_genr_str, res.abcd, 8)
  perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr_str,brain_genr_str,nperm=999, 
                                                    res.abcd,abs(cor.abcdTogenr))
  
  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr_total$pval.perm))
})

saveRDS(str_train_test_wisc_abcd, "str_train_test_wisc_abcd.rds")
setwd("/Users/estella/Desktop/Structural_cca/")
str_train_test_wisc_abcd <- readRDS("str_train_test_wisc_abcd.rds")
i <- 7
structure_resample <- str_train[[i]]
wisc <- pearson_abcd %>% filter(idc %in% id_train[[i]]) %>% mutate(wisc=as.numeric(pea_wiscv_tss))
wisc <- wisc[order(wisc$idc),]

cbcl_train <- structure_resample$cbcl_train
cbcl_train$wisc <- wisc$wisc
cbcl_wisc_train <- cbcl_train
print(nrow(cbcl_wisc_train))
id_na <- which(is.na(cbcl_wisc_train$wisc))
cbcl_wisc_train <- cbcl_wisc_train[-id_na,]
print(nrow(cbcl_wisc_train))
brain_train <- structure_resample$brain_train
brain_train <- brain_train[-id_na,]
print(nrow(brain_train))


structure_resample_test <- str_test[[i]]
wisc0 <- pearson_abcd %>% filter(idc %in% id_test[[i]]) %>% mutate(wisc=as.numeric(pea_wiscv_tss))
wisc0 <- wisc0[order(wisc0$idc),]

cbcl_test <- structure_resample_test$cbcl_test
cbcl_test$wisc <- wisc0$wisc
cbcl_wisc_test <- cbcl_test
print(nrow(cbcl_wisc_test))
id_na <- which(is.na(cbcl_wisc_test$wisc))
cbcl_wisc_test <- cbcl_wisc_test[-id_na,]
print(nrow(cbcl_wisc_test))
brain_test <- structure_resample_test$brain_test
brain_test <- brain_test[-id_na,]
print(nrow(brain_test))

brain_pen <- grid.abcd.str.wisc[[i]][[1]][[1]]
cbcl_pen <- grid.abcd.str.wisc[[i]][[1]][[2]]

res.abcd <- CCA(x=brain_train, z=cbcl_wisc_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                typex="standard", typez="standard",niter = 20, K=8)

# visualize the loadings:
cbcl_loading <- abs(res.abcd$v)
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression","IQ")

corrplot(t(cbcl_loading)[1:3,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 2, number.cex=1)



########################################################
############## Calculate canonical scores ##############
########################################################
CV_abcd_total <- scale(brain_whole) %*% brain_mean
cv1_abcd <- CV_abcd_total[,1]
cv2_abcd <- CV_abcd_total[,2]
cv3_abcd <- CV_abcd_total[,3]

df_cvs <- data.frame(participant_id=id, cv1=cv1_abcd, cv2=cv2_abcd, cv3=cv3_abcd)


########################################################
############## Merge the data ##########################
########################################################

all0 <- merge(df_cvs, pearson_abcd, by="participant_id")
dim(all0)
all1 <- merge(all0, wisc_abcd1, by="participant_id")
dim(all1)
all2 <- merge(all1, allinfo[,c("participant_id","sex","site","race_ethnicity","age","parental_education")], by="participant_id")
colSums(is.na(all2))
all_noNA <- all2
dim(all_noNA)
names(all_noNA)[5:8] <- c("matrix","fluid","cryst","total")
str(all_noNA)
all_noNA$matrix <- as.numeric(as.character(all_noNA$matrix))
all_noNA$fluid <- as.numeric(as.character(all_noNA$fluid))
all_noNA$cryst <- as.numeric(as.character(all_noNA$cryst))
all_noNA$total <- as.numeric(as.character(all_noNA$total))

all_noNA <- all_noNA[rowSums(is.na(all_noNA[,5:8])) == 0, ]
dim(all_noNA)
all_noNA <- all_noNA[all_noNA$cryst != max(all_noNA$cryst), ]# there's an outlier in crystal
dim(all_noNA)
########################################################
############## linear regressions ######################
########################################################
###### CV1
names(all_noNA)
matrix <- lm(scale(matrix) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(matrix)
confint(matrix)
fluid <- lm(scale(fluid) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(fluid)
confint(fluid)
cryst <- lm(scale(cryst) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(cryst)
confint(cryst)
total <- lm(scale(total) ~ scale(cv1)+site+sex+race_ethnicity+age+parental_education, all_noNA)
summary(total)
confint(total)


mat <- cor(df, use="complete.obs")
corrplot(mat)

summary(all_noNA)

ggplot(all_noNA,aes(cv1, cryst)) +
  geom_point() +
  geom_smooth(method='lm') 


p.adjust(c(0.10762,0.22264,0.220253,0.0262,
           0.00341,0.012060,0.000997,0.000183,
           0.0148,0.000105,0.000647,0.000001))


cormat <- readRDS("cormat.rds")

corrplot(cormat[1:6, -c(1:7)],method="color", 
         addCoef.col = "black", tl.srt =45, 
         tl.col = "black", tl.cex = 0.9, number.cex=0.8,
         col.lim = c(-0.2,0.2),
         is.corr = FALSE,
         col = colorRampPalette(c("blue", "white", "yellow"))(100))























########################################################
############## make the process automatic  #############
########################################################
structure_train_test <- function(brain_train,cbcl_train,brain_test,cbcl_test,brain_pen,cbcl_pen) {
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd, 8)
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=999,res.abcd,abcd.test)
  
  cor.abcdTogenr <- test_project_weights(brain_genr, cbcl_genr, res.abcd, 8)
  perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr,brain_genr,nperm=999, 
                                                    res.abcd,abs(cor.abcdTogenr))
  return(list(abcd.train=res.abcd,abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
         res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr_total$pval.perm))
}


########## results without re-ordering 
split_10_results <- lapply(1:10, function(i) {
  structure_resample <- str_resample[[i]]
  brain_train <- structure_resample$brain_train
  brain_test <- structure_resample$brain_test
  cbcl_train <- structure_resample$cbcl_train
  cbcl_test <- structure_resample$cbcl_test
          
  brain_pen <- grid.abcd[[i]][[1]][1]
  cbcl_pen <- grid.abcd[[i]][[1]][2]
        
  split <- structure_train_test(brain_train,cbcl_train,brain_test,cbcl_test,brain_pen,cbcl_pen)
    })

setwd("/Users/estella/Desktop/Structural_cca/")
saveRDS(split_10_results, "split_10_results.rds")

# plot the canonical loadings

cbcl_loading <- abs(split_10_results[[2]]$abcd.train$v)
rownames(cbcl_loading)  <- c("anxious","withdrawn","somatic","social","thought","attention","rule_breaking","aggression")

corrplot(t(cbcl_loading)[1:3,], method="color", 
         addCoef.col = "black", tl.srt =45, tl.col = "black", tl.cex = 2, number.cex=1)



##############################################################################
############ use split 3 to further explore the canonical loadings ###########
##############################################################################

strbrain_loading <- abs(split_10_results[[2]]$abcd.train$u)
strbrain_loading <- strbrain_loading[, 1:3]
rownames(strbrain_loading) <- colnames(brain_abcd)[-1]
colnames(strbrain_loading) <- rownames(c("CV1","CV2","CV3"))

str_brain_loadings <- strbrain_loading[rowSums(strbrain_loading) != 0, ]
corrplot(t(str_brain_loadings), method="color", 
         addCoef.col = "black", tl.srt =45, 
         tl.col = "black", tl.cex = 1.5, number.cex=1)

###########################################################################
###### bootstrap of the best split: stability check of cbcl and the brain
###########################################################################
set.seed(519)
######## bootstrap of the CVs and correlations
cl <- makePSOCKcluster(8)
registerDoParallel(cl)
n_boot <- 1000
boot_abcd <- foreach(1:n_boot, .packages="PMA") %dopar% {
  # boot samples
  N <- nrow(cbcl_train)
  idx <- sample.int(N, N, replace = TRUE) 
  # index the boot samples
  brain_train_boot <- brain_train[idx, ]
  cbcl_train_boot <- cbcl_train[idx, ]
  # run cca
  abcd_boot <- CCA(x=brain_train_boot, z=cbcl_train_boot, penaltyx = 0.6, penaltyz = 0.6, 
                   typex="standard", typez="standard",niter = 20, K=8)
  # extract loadings
  cbcl_loading <- abcd_boot$v
  brain_loading <- abcd_boot$u
  # test on the test set
  cor.abcd.test_boot <- test_project_weights(brain_test, cbcl_test, abcd_boot, 8)
  
  list(cbcl_loading=cbcl_loading, brain_loading=brain_loading, 
       cors_train=abcd_boot$cors, cors_test=cor.abcd.test_boot)
}
stopCluster(cl)

###### reorder the CVs of cbcl and brain
load_std <- res.abcd$u
cv_reorder_boot_abcd  <- lapply(seq_along(boot_abcd), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:4) {
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:4, function(z) {cor(load_std[, i], boot_abcd[[j]]$brain_loading[, z])})
    if(max(abs(cor_cv)) < 0.5){ # sometimes the correlation is too low
      idx  <- i
    } else {
      idx <- which(abs(cor_cv) == max(abs(cor_cv)))
    }
    
    if(!idx %in% idx_v) {
      idx_v <- c(idx_v, idx) # sometimes it will overlap
    } else {
      idx_v <- c(idx_v, i) 
    }
  }
  
  idx_v[duplicated(idx_v)] <- c(1:4)[!c(1:4) %in% idx_v]
  list(idx_reorder=idx_v)
  
})



# training set
boot_cor_train <- lapply(boot_abcd, function(x) {x$cors_train[1:4]})
boot_cor_train <- do.call(rbind,boot_cor_train)
boot_cor_train_reorder <- lapply(1:n_boot, function(i) {boot_cor_train[i, cv_reorder_boot_abcd[[i]]$idx_reorder]})
boot_cor_train_reorder <- do.call(rbind, boot_cor_train_reorder)
colnames(boot_cor_train_reorder) <- paste0("CV",1:4)
meantrain_boot <- colMeans(boot_cor_train)
sdtrain_boot <- apply(boot_cor_train, 2, sd)

# test set
boot_cor_test <- lapply(boot_abcd, function(x) {x$cors_test[1:4]})
boot_cor_test <- do.call(rbind,boot_cor_test)
boot_cor_test_reorder <- lapply(1:n_boot, function(i) {boot_cor_test[i, cv_reorder_boot_abcd[[i]]$idx_reorder]})
boot_cor_test_reorder <- do.call(rbind,boot_cor_test_reorder)
colnames(boot_cor_test_reorder) <- paste0("CV",1:4)
meantest_boot <- colMeans(boot_cor_test)
sdtest_boot <- apply(boot_cor_test, 2, sd)


fig_boot_syn <- data.frame(Train_Test = rep(c("Training set","Test set"),each=4), meancor = c(meantrain_boot, meantest_boot),
                           sdcor = c(sdtrain_boot, sdtest_boot), CV=rep(paste0("CV",1:4),2))


# violin plot with error bars
#########################
#train
df1 <- as.data.frame(boot_cor_train_reorder)
names(df1) <- paste0("CV", 1:4)
df1$traintest <- "Train"
#test
df2 <- as.data.frame(boot_cor_test_reorder)
names(df2) <- paste0("CV", 1:4)
df2$traintest <- "Test"

# combine train and test
df3 <- rbind(df1,df2)
df3 <- df3[, c(1:3,5)]
# reshape the data frame
df_long <- gather(df3, CV, cors, CV1:CV3, factor_key=TRUE)
df_long$cors <- abs(df_long$cors)
range(df_long[df_long$traintest == "Test", "cors"])


summary(df_long)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

p3 <- ggplot(df_long, aes(x = CV, y = abs(cors), fill = traintest)) +
  geom_point(aes(col=traintest), position=position_jitter(height=.05, width=.25),alpha = .3)+ 
  introdataviz::geom_split_violin(alpha = .8, trim = FALSE) +
  stat_summary(fun.data = data_summary,geom="pointrange", 
               position = position_dodge(.3)) + theme_bw() + 
  scale_x_discrete(name = "Canonical Variates", labels = paste0("CV", 1:3)) +
  scale_y_continuous(name = "Canonical correlations") +
  scale_fill_manual(name="Training Test Set", values = c("#7D3C98", "#52BE80")) + 
  scale_color_manual(name="Training Test Set", values = c("#BB8FCE", "#7DCEA0")) + 
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.title = element_text(size=15,face="bold"),
        legend.text = element_text(size=15),
        legend.position=c(0.25, 0.2)) 
p3

ggsave("train_test_boot_violin_reorder.pdf")


######################################################
####### visualization of the CIs of cbcl loadings ####
######################################################
cbclnames <- c("Anxious","Withdrawn","Somatic","Social","Thought","Attention","Rule Breaking","Aggression")
# training set
n_boot <- 1000
boot_cbcl <- lapply(1:n_boot, function(i) {
  boot_abcd[[i]]$cbcl_loading[,cv_reorder_boot_abcd[[i]]$idx_reorder]})

cbcl_sign <- sign(res.abcd$v)
cbcl_sign <- c(1,-1,1)


# cv1: attention, cv2: rule breaking, cv3:withdrawn, all +1

df_cbcl <- lapply(1:3, function(i) {
  cbcl_cv <- do.call(cbind, lapply(boot_cbcl, function(x){x[, i]}))
  cbcl_cv <- t(abs(cbcl_cv) * cbcl_sign[i])
  cbcl_cv <- as.data.frame(cbcl_cv)
  names(cbcl_cv) <- cbclnames
  cbcl_cv$CV <- paste0("CV", i)
  return(cbcl_cv)})

df_cbcl <- do.call(rbind,df_cbcl)
df_long_cbcl <- gather(df_cbcl, cbclsyndromes, loadings, cbclnames, factor_key = TRUE)
df_cbcl1 <- df_cbcl %>% filter(CV == "CV1")


############ combine three CVs together
# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND treatment :
pdf(file="cbcl_loadings_boot.pdf", width=13, height=5)
myplot <- boxplot(loadings ~ CV*cbclsyndromes, data=df_long_cbcl, 
                  boxwex=0.5, ylab="Canonical Loadings", outline=FALSE,
                  col=c("#FFDF6F","#B6DEE7","#D65353"), 
                  border=c("#F39C12","#2B78A2","#C0392B"),
                  xaxt="n", alpha=0.7,ylim = c(-1, 1))
# Add the grey vertical lines
for(i in seq(0.5 , 24, 3)){ 
  abline(v=i,lty=1, col="grey")
}
# To add the label of x axis
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 3)]
axis(1, font=2,
     at = seq(2, 24, 3), 
     labels = my_names , 
     tick=FALSE,lty = 2,lwd = 0.5)
# Add a legend
legend(22,-0.2, legend = c("CV1", "CV2", "CV3"), 
       col=c("#FFDF6F","#B6DEE7","#D65353"),
       border=c("#F39C12","#2B78A2","#C0392B"),
       pch = 15,  pt.cex = 2, cex = 0.9,  horiz = F, inset = c(0.1, 0.1))
dev.off()



######################################################
####### visualization of the CIs of brain loadings ####
######################################################
# training set
n_boot <- 1000

boot_brain <- lapply(1:n_boot, function(i) {
  boot_abcd[[i]]$brain_loading[,cv_reorder_boot_abcd[[i]]$idx_reorder]})
str_brain_loadings <- res.abcd$u
brain_sign <- sign(res.abcd$u)
max1 <- which(rank(-abs(res.abcd$u[,1])) == 1)
sign1 <- sign(res.abcd$u[max1,1])
max2 <- which(rank(-abs(res.abcd$u[,2])) == 1)
sign2 <- sign(res.abcd$u[max2,2])
max3 <- which(rank(-abs(res.abcd$u[,3])) == 1)
sign3 <- sign(res.abcd$u[max3,3])

sign <- c(sign1,sign2,sign3)

df_brain <- lapply(1:3, function(i) {
  brain_cv <- do.call(cbind, lapply(boot_brain, function(x){x[, i]}))
  brain_cv <- t((brain_cv))
  brain_cv <- as.data.frame(brain_cv)
  # adjust the signs 
  cvmean <- colMeans(abs(brain_cv))
  idx <- which(rank(-cvmean) == 1) 
  
  for (j in 1:nrow(brain_cv)){
    if(sign(brain_cv[j, idx]) != sign[i]){
      brain_cv[j,] <- -brain_cv[j, ]
    }
  }
   
  brain_cv$CV <- paste0("CV", i)
  return(brain_cv)})


df_brain <- do.call(rbind,df_brain)
names(df_brain)

colMeans(abs(df_brain[, -ncol(df_brain)]))
# inferiorparietal_vol, 7
# "posteriorcingulate_vol" 
# "precentral_vol" 
#  "supramarginal_vol"

# "caudalanteriorcingulate_surfarea"
# "postcentral_vol" 
# "posteriorcingulate_vol" 
# "rostralmiddlefrontal_vol" 
# "superiorfrontal_vol"
# "supramarginal_vol"


brain_check <- abs(res.abcd$u[,1:3])
quantile(brain_check,0.90)
brain_check[brain_check < 0.21] <- 0
rownames(brain_check) <- colnames(brain_train)
brain_check <- brain_check[rowSums(brain_check) != 0, ]
brain_check

###################################
########### plot the whole range
###################################
df1 <- df_brain[df_brain$CV == "CV1", 1:109]
df1 <- gather(df1, brain, loadings, bankssts_vol:Pallidum_vol, factor_key = TRUE)
dim(df1)
df1$measures <- c(rep(c("cortical volumes", "cortical surface areas", "cortical thickness"), 
                      each=34*1000), 
                  rep("subcortical volumes", 7*1000))

ggplot(df1, aes(x=brain, y = loadings)) + 
  geom_point(size = 0.8, aes(colour = measures)) + 
  scale_color_manual(values = c("#92C4EE", "#96D2CF", "#c5a0c5", "#B8DE9A")) + 
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face="bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15), 
        legend.position="bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))

setwd("/Users/estella/Desktop/Structural_cca/Figures")
ggsave("CV1_whole_brain.pdf", width = 15, height = 5)





# select the important brain structures
str_select <- c("transversetemporal_vol","supramarginal_vol","precentral_vol",
                "inferiorparietal_surfarea", "postcentral_surfarea",
                "superiorfrontal_surfarea","Amygdala_vol")

length(str)
df_brain1 <- df_brain[, names(df_brain) %in% c(str_select,"CV")]
#brain_structures <- str_select 
df_long_brain <- gather(df_brain1,brain_structures,loadings,str_select,factor_key = TRUE)


############ combine three CVs together
# Then I make the boxplot, asking to use the 2 factors : variety (in the good order) AND treatment :
pdf(file="brain_loadings_boot.pdf", width=13, height=5)

myplot <- boxplot(loadings ~ CV*brain_structures, data=df_long_brain, 
                  boxwex=0.5, ylab="Canonical Loadings", outline=FALSE,
                  col=c("#FFDF6F","#B6DEE7","#D65353"),
                  border=c("#F39C12","#2B78A2","#C0392B"),
                  xaxt="n", alpha=0.7,ylim = c(-1, 1))
# Add the grey vertical lines
for(i in seq(0.5, 27, 3)){ 
  abline(v=i,lty=1, col="grey")
}
# To add the label of x axis
my_names <- sapply(strsplit(myplot$names, '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 3)]
axis(1, at=seq(2,27,3),labels=FALSE)
text(x = seq(2,22,3), font=2,
     ## Move the labels down by 0.45.
     y = par("usr")[3]-0.1,
     labels = my_names,
     xpd = NA,
     srt = 25,
     cex = 0.9)
# Add a legend
legend(19.5,-0.4, legend = c("CV1", "CV2", "CV3"), 
       col=c("#FFDF6F","#B6DEE7","#D65353"),
       border=c("#F39C12","#2B78A2","#C0392B"),
       pch = 15,  pt.cex = 2, cex = 0.9,  horiz = F, inset = c(0.1, 0.1))

dev.off()














