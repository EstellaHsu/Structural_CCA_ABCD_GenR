#####################################################################################
################ bootstrap: stability check of cbcl and the brain ###################
#####################################################################################


# choose one train-test splits as the reference split (it can be done in any split)
brain_train <- str_train[[n]]$brain_train # n could be any split number 
cbcl_train <- str_train[[n]]$cbcl_train
brain_test <- str_test[[n]]$brain_test
cbcl_test <- str_test[[n]]$cbcl_test

res.abcd <- CCA(x=brain_train, z=cbcl_train, penaltyx = 0.6, penaltyz = 0.6, typex="standard", typez="standard",niter = 20, K=8)
# the penalty parameters are the ones you got in grid search for the specific split you chose

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
                   # the penalty parameters are the ones you got in grid search for the specific split you chose
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
load_std <- res.abcd$u # based on the brain loadings
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

df_cbcl <- lapply(1:3, function(i) {
  cbcl_cv <- do.call(cbind, lapply(boot_cbcl, function(x){x[, i]}))
  cbcl_cv <- t(abs(cbcl_cv) * cbcl_sign[i])
  cbcl_cv <- as.data.frame(cbcl_cv)
  names(cbcl_cv) <- cbclnames
  cbcl_cv$CV <- paste0("CV", i)
  return(cbcl_cv)})

df_cbcl <- do.call(rbind,df_cbcl)
df_long_cbcl <- gather(df_cbcl, cbclsyndromes, loadings, cbclnames, factor_key = TRUE)

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

ggsave("CV1_whole_brain.pdf", width = 15, height = 5)



# select the important brain structures
str_select <- c("transversetemporal_vol","supramarginal_vol","precentral_vol",
                "inferiorparietal_surfarea", "postcentral_surfarea",
                "superiorfrontal_surfarea","Amygdala_vol")

length(str)
df_brain1 <- df_brain[, names(df_brain) %in% c(str_select,"CV")]
brain_structures <- str_select 
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




