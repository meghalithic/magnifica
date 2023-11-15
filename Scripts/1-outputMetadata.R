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

output.fossil <- read.csv("./Data/output_4Aug2023_done.csv", header = TRUE)
nrow(output.fossil) #6443

output.modern <- read.csv("./Data/output_modern_trim.csv",
                          header = TRUE)
nrow(output.modern) #922
## remove images that are not magnifica
colnames(output.modern)
head(output.modern$id)
#images are 1200, 1203, 1205, 1216, 1220
rm.img <- c("bryozoa_lab_images/mab-modern-jpg/1200_CC_1_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1200_CC_2_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1203_CC_1_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1203_CC_2_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1203_CC_3_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1205_CC_1_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1205_CC_2_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1205_CC_3_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1205_CC_4_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1216_CC_1_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1216_CC_2_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1216_CC_3_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1220_CC_1_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1220_CC_2_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1220_CC_3_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1220_CC_4_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1220_CC_5_15v_x30_BSE.jpg",
            "bryozoa_lab_images/mab-modern-jpg/1220_CC_6_15v_x30_BSE.jpg")
output.modern <- output.modern[!(output.modern$id %in% rm.img),]
nrow(output.modern) #613

bryo.meta <- read.csv("./Data/image_merge_txt_usingfileName_DONE_17Apr2023.csv",
                      header = TRUE,
                      sep = ";")

#### EXPLORE DATA ----

nrow(bryo.meta) #1890

output.modern$image <- gsub("bryozoa_lab_images/mab-modern-jpg/",
                             "",
                             output.modern$id)


output.fossil$image <- gsub("bryozoa_lab_images/",
                             "",
                             output.fossil$id)

output <- as.data.frame(rbind(output.fossil, output.modern))


output$image <- gsub(".jpg",
                     "",
                     output$image)

nrow(output) #7056
colnames(output)
output$id[1]

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
                     by.x = "image", by.y = "imageName", #don't use image for bryo.meta.trim
                     all.x = TRUE, all.y = FALSE)

nrow(meta.images) #7056 #zooids
length(unique(meta.images$image)) #1455 images

unique(meta.images$formation)
meta.images$image[is.na(meta.images$formation)]
meta.images$formation[is.na(meta.images$formation)] <- "modern"
unique(meta.images$formation)

sp.nr <- meta.images$image[is.na(meta.images$newSpecimenNR)]
sp.nr.sp <- str_split(sp.nr, "_")
sp.nr.new <- c()
for(i in 1:length(sp.nr.sp)){
    sp.nr.new[i] <- paste0(sp.nr.sp[[i]][1], sp.nr.sp[[i]][2])
}
meta.images$newSpecimenNR[is.na(meta.images$newSpecimenNR)] <- sp.nr.new

length(unique(meta.images$newSpecimenNR)) #751 colonies

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
          "./Data/meta.images_15Nov2023.csv", #includes modern data now
          row.names = FALSE)

