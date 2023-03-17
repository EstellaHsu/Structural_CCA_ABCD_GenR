########### data preparation ##############
###########################################

library(data.table)
library(dplyr)
library(stringr)

###########################
#structural data QC
###########################
setwd("V:/medewerkers/051950 Xu, B/PhD projects/Structural_CCA/abcd")
qc <- readRDS("all_final_abcd_tbv.rds")
dim(qc)

###########################
#structural data processing
###########################

cortical <- readRDS("abcd_freesurfer_20191219_aparc_stats.rds")
cortical <- cortical[cortical$idc %in% qc$idc, ]
subcortical <- readRDS("abcd_freesurfer_20191219_aseg_stats.rds")
subcortical <- subcortical[subcortical$idc %in% qc$idc, ]
names(cortical)
names(subcortical)
dim(cortical)
dim(subcortical)
identical(cortical$idc,subcortical$idc)

# function of average left and right cortical regions
ave_lhrh <- function(measures) {
  cor <- cortical %>% select(contains(measures)) %>% as.data.frame
  names(cor) <- names(cor) %>% gsub(c("lh_|rh_|_abcd"),"",.)
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


######## subcortical volumes

# subcortical areas of interest
sub_area <- c("Hippocampus","Amygdala","Accumbens","Caudate","Thalamus","Putamen","Pallidum")

sub1 <- lapply(sub_area, function(x) rowMeans(subcortical[, str_detect(names(subcortical),x)]))# the region names
sub2 <- do.call(cbind, sub1)
colnames(sub2) <- paste0(sub_area, "_vol")
subcor_vol <- as.data.frame(sub2)


###########################
#data merging
###########################
brain <- cbind(idc=cortical$idc, cor_vol, cor_surfarea, cor_thick, subcor_vol)
str(brain)
dim(brain)

df_all <- merge(qc, brain, by="idc")
colSums(is.na(df_all))

saveRDS(df_all,"all_final_abcd_tbv.rds")

names(df_all)

str(df_all)
###########################
#standardize the outcomes
###########################
df_all[23:131] <- lapply(df_all[23:131], scale)
all_abcd_zscores <- df_all
names(all_abcd_zscores)


df_all[22:130] <- lapply(df_all[22:130], scale)
all_abcd_zscores <- df_all
saveRDS(all_abcd_zscores, "all_final_abcd_std_13.rds")

