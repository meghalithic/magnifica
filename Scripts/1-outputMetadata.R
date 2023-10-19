## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) compares image list in output from Steginator to bryozoa metadata file

## the output of this file is meta.images_date.csv

#### LOAD DATA ----

source("./Scripts/0-env.R")

#output.oldest <- read.csv("./Data/output_May2023_Stegniator.csv", header = TRUE)
#nrow(output.oldest) #19346

#output.old <- read.csv("./Data/output_21Jun2023_Stegniator.csv", header = TRUE)
#nrow(output.old) #15783 #didn't ignore zooids, just filtered out images

#output <- read.csv("./Data/output_31jul2023_Steginator.csv", header = TRUE)
#nrow(output) #7202 

output <- read.csv("./Data/output_4Aug2023_done.csv", header = TRUE)
nrow(output) #6443

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
#882_1405_433_617
#598_1024_434_666
#880_646_517_625

output$new.id <- paste0(output$box_id, "_", output$image)
xx <- output[duplicated(output$new.id),]
nrow(xx) #0! yay!

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

####### COMPARE TO OLD SAMPLING NUMBERS ------

output.oldest$fileName <- gsub("/media/voje-lab/bryozoa/imaging/Stegino_images/combined/", "", output.oldest$id)

output.oldest$image <- gsub(".tif", "", output.oldest$fileName)

zz.old <- output.oldest[duplicated(output.oldest$box_id),]
nrow(zz.old) #3
zz.old
#82_1405_433_617
#598_1024_434_666
#880_646_517_625

output.oldest$new.id <- paste0(output.oldest$box_id, "_", output.oldest$image)
xx.old <- output.oldest[duplicated(output.oldest$new.id),]
nrow(xx.old) #0! yay!

meta.images.old <- merge(output.oldest, bryo.meta.trim,
                     by.x = "image", by.y = "imageName",
                     all.x = TRUE, all.y = FALSE)

nrow(meta.images.old) #19346
length(unique(meta.images.old$newSpecimenNR)) #904 colonies
length(unique(meta.images.old$fileName)) #1880 iimages

#### WRITE OUT DATASET ----
## CHANGE DATE EXT EVERYTIME
write.csv(meta.images,
          "./Data/meta.images_29Sept2023.csv",
          row.names = FALSE)

