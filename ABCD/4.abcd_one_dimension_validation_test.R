####################################################################################################
####################### ABCD post hoc analyses: one dimensional structure ##########################
####################################################################################################


########## load all the packages
packages <- c('psych','doParallel','permute','reshape2','PMA','caret','corrplot','groupdata2',
              'data.table','readxl','dplyr','stringr','base','Amelia','parallel','matrixStats',
              'gridExtra','grid','ggplot2','lattice','tidyverse','hrbrthemes','viridis','forcats',
              'sjPlot','lme4','factoextra', 'fmsb', 'NbClust','stats','dendextend','cluster',
              'fpc', 'e1071', 'plotly', 'MASS', 'ggradar2','tidyr','introdatavizyes','Hmisc',
              'pheatmap','RColorBrewer','circlize','CCA','CCP','candisc')


lapply(packages, require, character.only = TRUE)


####################################################################
############## traditional CCA with orthogonality ##################
####################################################################

##### We did the traditional CCA in one train-test split (it can be done in any split)
n <- sample(1:10, 1)
brain_train <- str_train[[n]]$brain_train
cbcl_train <- str_train[[n]]$cbcl_train
brain_test <- str_test[[n]]$brain_test
cbcl_test <- str_test[[n]]$cbcl_test

##### fit the traditional CCA model
abcd_train_cca <- candisc::cancor(brain_train, cbcl_train)


##### calculate the canonical correlations in the ABCD test set
std_brain <- scale(brain_test) %*% abcd_train_cca$coef$X
std_cbcl <- scale(cbcl_test) %*% abcd_train_cca$coef$Y

cor.res.test.abcd <- diag(cor(std_brain, std_cbcl))

##### test the significance using permutation tests
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


##### test the correlations in generation R
std_brain_genr <- scale(brain_genr_str) %*% abcd_train_cca$coef$X
std_cbcl_genr <- scale(cbcl_genr_str) %*% abcd_train_cca$coef$Y
cor.res.genr <- diag(cor(std_brain, std_cbcl))
cor.res.genr



##################################################################################
############## using item-level scores of CBCL in the CCA model ##################
##################################################################################
# NOTE: the data extraction procedure is the same with the CBCL syndrome scale scores, the only difference is the behavioral measures

########################################################
############## 1. Read the data ########################
########################################################
all_abcd_item <- readRDS("train_test_split_item.rds")

########################################################
############## 2. split in training and test ###########
########################################################

# extract the data first 
str_train_item <- lapply(1:10, function(i) {temp <- all_abcd_item[[i]]
list(brain_train=temp$brain_train, 
     cbcl_train=temp$cbcl_train)}) 

str_test_item <- lapply(1:10, function(i) {temp <- all_abcd_item[[i]]
list(brain_test=temp$brain_test, 
     cbcl_test=temp$cbcl_test)})

##### IMPORTANT! for the item level scores, some are with 0 variance, which is not allowed in the SCCA model, so we need to remove them
zerovar_test <- unique(unlist(lapply(str_test_item, function(x) {
  a <- apply(x$cbcl_test, 2, var)
  names(a)[which(a == 0)]})))

zerovar_train <- unique(unlist(lapply(str_train_item, function(x) {
  a <- apply(x$cbcl_train, 2, var)
  names(a)[which(a == 0)]})))

##### The items have 0 variance in all ABCD test sets are overlapped with ABCD training sets


##### Extract the data we need 
str_train_item <- lapply(1:10, function(i) {temp <- all_abcd_item[[i]]
list(brain_train=temp$brain_train[,c(1:68,103:109)], # index of brain measures, no cortical thickness
     cbcl_train=temp$cbcl_train[, !names(temp$cbcl_train) %in% zerovar_test])}) # remove the 0 variance items

str_test_item <- lapply(1:10, function(i) {temp <- all_abcd_item[[i]]
list(brain_test=temp$brain_test[,c(1:68,103:109)], # index of brain measures, no cortical thickness
     cbcl_test=temp$cbcl_test[, !names(temp$cbcl_test) %in% zerovar_test])}) # remove the 0 variance items


########################################################
############## 3. Fit the sCCA model  ##################
########################################################
###### The grid search of penalty parameters
grid.abcd.str.item <- lapply(1:10, function(i) {
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



###### Fit the model in ABCD and Generation R
str_train_test_abcd_item <- lapply(1:10, function(i) {
  brain_train <- str_train[[i]]$brain_train
  cbcl_train <- str_train[[i]]$cbcl_train
  brain_test <- str_test[[i]]$brain_test
  cbcl_test <- str_test[[i]]$cbcl_test
  brain_pen <- grid.abcd.str.item[[i]][1]
  cbcl_pen <- grid.abcd.str.item[[i]][2]
  
  res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = brain_pen, penaltyz = cbcl_pen, 
                  typex="standard", typez="standard",niter = 20, K=8)
  perm_abcd_train <- permutation_test(cbcl_train, brain_train, nperm=1999, cbcl_pen, brain_pen, 8, res.abcd$cors)
  
  abcd.test <- test_project_weights(brain_test,cbcl_test,res.abcd, 8)
  perm_abcd_test <- permutation_test_testset(cbcl_test,brain_test,nperm=1999,res.abcd,abcd.test)
  
  cor.abcdTogenr <- test_project_weights(brain_genr_str, cbcl_genr_str, res.abcd, 8)
  perm_abcdTogenr_total <- permutation_test_testset(cbcl_genr_str,brain_genr_str,nperm=1999, 
                                                    res.abcd,abs(cor.abcdTogenr))
  
  return(list(abcd.train=res.abcd,abcd.train.perm=perm_abcd_train$pval.perm,
              abcd.test=abcd.test,abcd.test.perm=perm_abcd_test$pval.perm,
              res.genr=cor.abcdTogenr,genr.perm=perm_abcdTogenr_total$pval.perm))
})


saveRDS(str_train_test_abcd_item,"str_train_test_abcd_item.rds")



##########################################################################################
############## below are all about visualization of item-level results  ##################
##########################################################################################


########################################################
############## reorder and average across splits  ######
########################################################
# set the reference split (you can randomly choose one or select the one that you prefer the orders)
res.abcd.item <- CCA(x=brain_train, z=cbcl_train, penaltyx = 0.8, penaltyz = 0.4, typex="standard", typez="standard",
                niter = 20, K=8)


# the reference split
load_std <- res.abcd.item$v # reorder based on brain loadings
# we only care about the first 4 CVs based on the covariance explained plot. 
# As the CBCL loadings are high similar, if we have more CVs, the brain might be chaotic

new_cv_reorder  <- lapply(seq_along(str_train_test_abcd_item), function(j) {
  # go through all the splits
  idx_v <- c()
  for (i in 1:5) { # we specified 5 CVs here after inspecting results from several train-test splits, only the first 5 are relevant. 
    # calculate the maximum correlation of each CV and reorder
    cor_cv <- sapply(1:5, function(z) {cor(load_std[, i], str_train_test_abcd_item[[j]]$abcd.train$v[, z])})
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
  
  loadMat <- str_train_test_abcd_item[[j]]$abcd.train$v[, idx_v]
  return(list(loadMat=loadMat, idx_reorder=idx_v))
  
})

##### average loading of CBCL across 10 splits
cbcl <- lapply(new_cv_reorder, function(x) {x$loadMat})
cbcl_mean <- do.call(cbind, lapply(1:5, 
                                   function(i) {
                                     rowMeans(abs(do.call(cbind, 
                                                          lapply(cbcl, function(x) {x[, i]})
                                     )))}))

rownames(cbcl_mean)  <- items_name

##### average loading of brain structures across 10 splits

brain <- lapply(seq_along(new_cv_reorder), function(x) {
  loadMat <- str_train_test_abcd_item[[x]]$abcd.train$u[, new_cv_reorder[[x]]$idx_reorder]})

brain_mean <- do.call(cbind, lapply(1:5, 
                                    function(i) {
                                      rowMeans(abs(do.call(cbind, lapply(brain, function(x) {x[, i]})
                                      )))}))

rownames(brain_mean) <- colnames(brain_train)


########################################################
############## visualization of behavior and brain  ####
########################################################

###########################################
###### brain figures: cortical ############
###########################################

# only visualize the first 3 canonical variates
brain_fig <- brain_mean[,1:3]
# choose the 75% quantile
brain_fig <- round(brain_fig, 3)

rownames(brain_fig) <- colnames(brain_train)
colnames(brain_fig) <- paste0("CV",1:3)

####### Below is an example of creating the figure for surface areas (CV1)
####### Others can be created using the same code

# create a figure data frame of surface areas
fig_surf <- brain_fig[35:68, ] # index of surface areas
q_25 <- quantile(fig_surf) # choose the top 25% loading
fig_surf[fig_surf < q_25] <- 0
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


######## 2D figures
ggseg(.data=fig_surf[,c("region","CV1")], mapping=aes(fill=CV1),colour="grey") +
  labs(title="CV1 Surface areas", fill="loadings") +
  scale_fill_gradient(low="white",high="#FFDF6F",na.value="#F2F3F4")

ggsave("CV1_surfarea_item.pdf",width=10,height=4)


###########################################
###### brain figures: subcortical #########
###########################################
# create a figure data frame of subcortical volumes
fig_sub <- brain_fig[69:75, ] # index of subcortical volumes
q_25 <- quantile(fig_sub) # choose the top 25% loading
fig_sub[fig_sub < q_25] <- 0


####### Below is an example of creating the figure for CV1
cv_loadings <- fig_sub[,1]
subcor <- gsub("_vol","",rownames(fig_sub))
subcor <- c("hippocampus", "amygdala","accumbens area","caudate","thalamus proper","putamen","pallidum")
fig_str <- data.frame(region = subcor, loadings=fig_sub[,1])

fig_aseg <-  aseg_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  dplyr::select(label) %>% 
  filter(grepl("Hippocampus|Amygdala|Accumbens-area|Caudate|Thalamus-Proper|Putamen|Pallidum", 
               label)) %>% 
  mutate(loadings = rep(cv_loadings, 2)) 

p <- ggseg3d(.data = fig_aseg, atlas = aseg_3d, 
             colour = "loadings", text = "loadings", 
             na.alpha= .5,
             palette = c("white" = 0, "#FFDF6F" = .075, "#EFB603" = 0.15)) %>% 
  remove_axes() %>%
  add_glassbrain(colour = "white",
                 opacity = 0.2)
p

########################################
###### CBCL circular barplot ###########
########################################

cbcl_mean_syn  <- as.data.frame(cbcl_mean)
cbcl_mean_syn <- cbcl_mean[items_syndrome, 1:3]

quantile(cbcl_mean[,1:3], 0.90)
cbcl_mean_syn[abs(cbcl_mean_syn) < 0.111279] <- 0
syndrome <- c(rep("attention", 10), rep("aggression", 18), rep("rule_breaking", 13),
              rep("anxious", 13), rep("withdrawn", 8), rep("somatic", 11),
              rep("social", 11), rep("thought", 15), rep("others",16))


fig_cbcl_item <- data.frame(items = items_syndrome, syndrome = syndrome, 
                            CV1=cbcl_mean_syn[, 1], CV2=cbcl_mean_syn[,2], CV3=cbcl_mean_syn[,3])

fig_cbcl_item <- fig_cbcl_item[rowSums(fig_cbcl_item[,c("CV1","CV2","CV3")]) != 0, ]

df_long_cbcl_item <- gather(fig_cbcl_item,canonical_variates,loadings,CV1:CV3,factor_key = TRUE)


data <- df_long_cbcl_item
# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 2
nObsType <- nlevels(as.factor(data$canonical_variate))
to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(data$syndrome))*nObsType, ncol(data)))
colnames(to_add) <- colnames(data)
to_add$syndrome <- rep(levels(as.factor(data$syndrome)), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(syndrome, items)
data$id <- rep(seq(1, nrow(data)/nObsType) , each=nObsType)

# Get the name and the y position of each label
label_data <- data %>% group_by(id, items) %>% summarise(tot=sum(loadings))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substracted 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(syndrome) %>% 
  summarise(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

p <- ggplot(data) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=loadings, fill=canonical_variates), 
           stat="identity", alpha=0.9) +
  scale_fill_manual(values = c("#FFDF6F","#B6DEE7","#D65353")) +
  
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.05, xend = start, yend = 0.05), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.10, xend = start, yend = 0.10), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.15, xend = start, yend = 0.15), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.20, xend = start, yend = 0.20), 
               colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(max(data$id),5), y = c(0, 0.05, 0.1, 0.15, 0.2), 
                    label = c("0", "0.05", "0.10", "0.15", "0.20") , color="grey",
                    size=3, angle=0, fontface="bold", hjust=1) +
  ylim(-0.2,0.5) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=tot + 0.01, label=items, hjust=hjust), color="black", 
            fontface="bold",alpha=0.6, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start-0.5, y = -0.01, xend = end+0.5, yend =-0.01), colour = "grey", 
               alpha=0.8, size=1.2, inherit.aes = FALSE)  
#geom_text(data=base_data, aes(x = title, y = -0.02, label=syndrome), hjust=1, 
#colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)
# Save at png
ggsave(p, file="cbcl_item.pdf", width=10, height=10)


