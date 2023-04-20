## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) reduces image list to those with 30 magnification

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)
require(ggplot2)
require(reshape2)
require(lmodel2)
require(tidyverse)
require(stringr)

#### LOAD DATA ----
bryo.images <- read.table("./Data/Steginoporella_magnifica_image_metadata_17Apr2023.csv",
                          header = TRUE,
                          sep = ";")
nrow(bryo.images) #1890

bryo.images.30 <- bryo.images[bryo.images$Magnification == 30,]
nrow(bryo.images.30) #1835

##### add in full file name #####
bryo.images.30$fileName.tif <- paste0(bryo.images.30$fileName, ".tif")

#### CREATE OLD FILE NAMES - DELETE LATER ####
recon.df <- read.table("./Data/image_merge_txt_usingfileName_DONE_17Apr2023.csv",
                       header = TRUE, sep = ";")

setdiff(recon.df$newFileName, bryo.images$fileName) #none
setdiff(bryo.images$fileName, recon.df$newFileName) #none

##### combine files to get path #####
recon.df.path <- recon.df[, c("newFileName", "path.tif")]
bryo.images.path <- merge(bryo.images.30, recon.df.path,
                          by.x = "fileName", by.y = "newFileName",
                          all.x = TRUE, all.y = FALSE)
nrow(bryo.images.path)
#duplicates to reconcile manually

#write.csv(bryo.images.path,
#          "./Data/filteredImages.csv",
#          row.names = FALSE)
