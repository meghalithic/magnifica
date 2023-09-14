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

output.oldest <- read.csv("./Data/output_May2023_Stegniator.csv", header = TRUE)
nrow(output.oldest) #19346

#output.old <- read.csv("./Data/output_21Jun2023_Stegniator.csv", header = TRUE)
#nrow(output.old) #15783 #didn't ignore zooids, just filtered out images

#output <- read.csv("./Data/output_31jul2023_Steginator.csv", header = TRUE)
#nrow(output) #7202 

output <- read.csv("./Data/output_4Aug2023_done.csv", header = TRUE)
nrow(output) #6443


#AP_images <- read.csv("./Data/images_from_AP.csv", header = TRUE)
bryo.meta <- read.csv("./Data/image_merge_txt_usingfileName_DONE_17Apr2023.csv",
                      header = TRUE,
                      sep = ";")

#### EXPLORE DATA ----

nrow(bryo.meta) #1890

nrow(output) #6443
colnames(output)
output$id[1]

output$fileName <- gsub("bryozoa_lab_images/",
                        "",
                        output$id)

output$image <- gsub(".jpg",
                     "",
                     output$fileName)

zz <- output[duplicated(output$box_id),]
nrow(zz) #3
zz
#82_1405_433_617
#598_1024_434_666
#880_646_517_625

output$new.id <- paste0(output$box_id, "_", output$image)
xx <- output[duplicated(output$new.id),]
nrow(xx) #0! yay!

output.oldest$fileName <- gsub("/media/voje-lab/bryozoa/imaging/Stegino_images/combined/",
                        "",
                        output.oldest$id)

output.oldest$image <- gsub(".tif",
                     "",
                     output.oldest$fileName)

zz.old <- output.oldest[duplicated(output.oldest$box_id),]
nrow(zz.old) #3
zz.old
#82_1405_433_617
#598_1024_434_666
#880_646_517_625

## DO WITH OLDER DATA SET SO THAT HAVE SOMETHING TO COMPARE
output.oldest$new.id <- paste0(output.oldest$box_id, "_", output.oldest$image)
xx.old <- output.oldest[duplicated(output.oldest$new.id),]
nrow(xx.old) #0! yay!


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

nrow(meta.images) #6443 #zooids
length(unique(meta.images$newSpecimenNR)) #732 colonies
length(unique(meta.images$fileName)) #1395 images

####### OLD SAMPLING NUMBERS ------

meta.images.old <- merge(output.oldest, bryo.meta.trim,
                     by.x = "image", by.y = "imageName",
                     all.x = TRUE, all.y = FALSE)

nrow(meta.images.old) #19346
length(unique(meta.images.old$newSpecimenNR)) #904 colonies
length(unique(meta.images.old$fileName)) #1880 iimages


## CHANGE DATE EXT EVERYTIME
write.csv(meta.images,
          "./Data/meta.images.8Sept2023.csv",
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
