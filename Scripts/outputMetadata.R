## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) compares image list in output from Steginator to bryozoa metadata file

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)
require(ggplot2)
require(reshape2)
require(lmodel2)
require(tidyverse)

#### LOAD DATA ----
## CHANGE DATA FILE NAME AS NEEDED
output <- read.csv("./Data/output_May2023_Stegniator.csv", header = TRUE)
#AP_images <- read.csv("./Data/images_from_AP.csv", header = TRUE)
bryo.meta <- read.csv("./Data/image_merge_txt_usingfileName_DONE_17Apr2023.csv",
                      header = TRUE,
                      sep = ";")

#### EXPLORE DATA ----

nrow(bryo.meta) #1890

nrow(output) #19346
colnames(output)
output$id[1]

output$fileName <- gsub("/media/voje-lab/bryozoa/imaging/Stegino_images/combined/",
                        "",
                        output$id)

output$image <- gsub(".tif",
                     "",
                     output$fileName)

zz <- output[duplicated(output$box_id),]
nrow(zz) #19

output$new.id <- paste0(output$box_id, "_", output$image)
xx <- output[duplicated(output$new.id),]
nrow(xx) #0! yay!

#output$specimenNR <- c()
#for(i in 1:nrow(output)){
#  output$specimenNR[i] <- paste0(str_split(output$image[i], fixed("_"))[[1]][1],
#                                 str_split(output$image[i], fixed("_"))[[1]][2])
#}

#### MERGE FILES ----

bryo.meta.trim <- bryo.meta %>%
  dplyr::select(newFileName, newImage, newSpecimenNR,
         image, fileName.tif, path.tif, specimenNR.tif,
         Enterer, NOTES, formation, Mag)

bryo.meta.trim$imageName <- gsub(".tif", "",
                                 bryo.meta.trim$fileName.tif)

meta.images <- merge(output, bryo.meta.trim,
                     by.x = "image", by.y = "imageName",
                     all.x = TRUE, all.y = FALSE)

nrow(meta.images) #19346


## CHANGE DATE EXT EVERYTIME
write.csv(meta.images,
          "./Data/meta.images.Jun2023.csv",
          row.names = FALSE)

#### OLD ----

#AP_meta <- read.table("./Data/AP_Steginoporella.csv", 
#                      header = TRUE,
#                      sep = ";")

##images from shared "JPG" folder
#nrow(AP_images) #1654
#imageName.parse_AP <- str_split(AP_images$imageName, fixed("_"))
#specimen.NR_AP <- c()
#for(i in 1:length(imageName.parse_AP)){
#  specimen.NR_AP[i] <- paste0(imageName.parse_AP[[i]][1], imageName.parse_AP[[i]][2])
#}
#AP_images$SPECIMEN.NR <- specimen.NR_AP
#nrow(AP_images[!duplicated(AP_images$SPECIMEN.NR),]) #779
#AP_images.trimmed <- AP_images[!duplicated(AP_images$SPECIMEN.NR),]

## overlap between A. Porto image list with bryo metadata:

#length(setdiff(AP_images.trimmed$SPECIMEN.NR, bryo.meta$SPECIMEN.NR)) #136
#length(setdiff(bryo.meta$SPECIMEN.NR, AP_images.trimmed$SPECIMEN.NR)) #134

## overlap between A. Porto image list with metadata provided by him:
length(setdiff(AP_images.trimmed$SPECIMEN.NR, AP_meta$specimen.nr.)) #137
length(setdiff(AP_meta$specimen.nr., AP_images.trimmed$SPECIMEN.NR)) #6

## overlap between bryo metadata and metadata provide by AP:
length(setdiff(bryo.meta$SPECIMEN.NR, AP_meta$specimen.nr.)) #0
##THIS IS GOOD! 

## overlap between A. Porto image list with output:
output.spID <- unique(output$SPECIMEN.NR)
jpg.ID <- unique(AP_images$SPECIMEN.NR)
setdiff(output.spID, jpg.ID) #"NA"
setdiff(jpg.ID, output.spID) #"3131.jpg" "492CC"    "494CC" not in output
#there is an image named 313_1.jpg; AP said to ignore
#492 and 494 look poor quality

##new match up csv between images and bryo meta
extra.images.AP <- setdiff(AP_images$SPECIMEN.NR, bryo.meta$SPECIMEN.NR)
extra.meta <- setdiff(bryo.meta$SPECIMEN.NR, AP_images$SPECIMEN.NR)

fromDataset <- c(rep("AP.images", length(extra.images.AP)), rep("metadata", length(extra.meta)))
imageDiff <- c(extra.images.AP, extra.meta)
AP.images.errors.df <- cbind(fromDataset, imageDiff)

##new match up csv between metadata from AP and bryo meta
extra.meta <- setdiff(bryo.meta$SPECIMEN.NR, AP_meta$specimen.nr.)
extra.meta.AP<- setdiff(AP_meta$specimen.nr., bryo.meta$SPECIMEN.NR)

fromDataset <- c(rep("AP.meta", length(extra.meta.AP)), rep("metadata", length(extra.meta)))
imageDiff <- c(extra.meta.AP, extra.meta)
meta.errors.df <- cbind(fromDataset, imageDiff)

#write.csv(AP.images.errors.df,
#          "./Results/AP.images.errors.df.csv",
#          row.names = FALSE)
## AP images are a subset of the metadata file, so it makes sense that there is a discrepancy
##

#write.csv(meta.errors.df,
#          "./Results/meta.errors.df.csv",
#          row.names = FALSE)
## AP metadata file is old, so makes sense there is a discrepancy
## the two images "not" in AP metadata that are found in the bryo metadata are actually there
