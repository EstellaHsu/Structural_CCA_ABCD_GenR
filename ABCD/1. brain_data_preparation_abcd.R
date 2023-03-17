############################################################
########### brain strcutural data preparation ##############
############################################################
library(data.table)
library(dplyr)
library(stringr)

###########################
#structural data QC
###########################
# first, read the ids that passed the QC, this is from previous step 0.
setwd("V:/medewerkers/***** Xu, B/PhD projects/Structural_CCA/abcd")
qc <- readRDS("all_final_abcd.rds")
dim(qc)

###########################
#structural data extraction
###########################

# extract the cortical data
cortical <- readRDS("abcd_freesurfer_20191219_aparc_stats.rds")
cortical <- cortical[cortical$idc %in% qc$idc, ]

# extract the subcortcial data
subcortical <- readRDS("abcd_freesurfer_20191219_aseg_stats.rds")
subcortical <- subcortical[subcortical$idc %in% qc$idc, ]
names(cortical)
names(subcortical)
dim(cortical)
dim(subcortical)
identical(cortical$idc,subcortical$idc)

# a function of averaging left and right cortical regions
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

# subcortical areas of interest: we selected 7 subcortical structures
sub_area <- c("Hippocampus","Amygdala","Accumbens","Caudate","Thalamus","Putamen","Pallidum")

# average the left and right hemisphere
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
# check whether there are NAs
colSums(is.na(df_all))
# No NAs

saveRDS(df_all,"all_final_abcd.rds")
