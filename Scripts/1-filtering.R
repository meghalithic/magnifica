
# Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) adds in formation info
# 2) filters images
# 3) makes sure that all datasets have formation and specimenNR

## The output of this code is images.filtered_8Dec2023.csv

#### ENVIRONMENT ----

source("./Scripts/0-env.R")

# # Manually set formations
# NKBS: 001-399
# NKLS: 400-599
# Tewkesbury: 600-699
# SHCSBSB: 700-799
# Tainui: 800-899
# Upper Kai-Iwi: 1000-1099
# Waipuru: 1100-1199
# Modern:1200-1299

#### LOAD DATA ----

### FILTERING TABLES
fossil.filter <- read.csv("Data/filteredImages.csv", sep = ";", header = TRUE)
bleed.filter <- read.csv("Data/BLEED.image.filter.csv", header = TRUE)
#keeping all wabo iv, so no filter file needed
#manually filtering modern, so no filter file needed

### FOSSIL
output.fossil <- read.csv("./Data/output_4Aug2023_done.csv", header = TRUE)
nrow(output.fossil) #6443
colnames(output.fossil)
head(output.fossil$id) #path to image

output.waboiv_1 <- read.csv("Data/output_waboiv.csv", header = TRUE)
output.waboiv_2 <- read.csv("Data/output.nkls.2.csv", header = TRUE)

nrow(output.waboiv_1) #370
colnames(output.waboiv_1)
head(output.waboiv_1$id) #path to image

nrow(output.waboiv_2) #224
colnames(output.waboiv_2)
head(output.waboiv_2$id) #path to image


#### COMBINE WABO IV
output.waboiv <- as.data.frame(rbind(output.waboiv_1, output.waboiv_2))
nrow(output.waboiv) #3594
colnames(output.waboiv)
head(output.waboiv$id) #path to image
## keep separate from output.fossil because id path is different

### MODERN

output.bleed <- read.csv("Data/output.bleed.csv", header = TRUE)

nrow(output.bleed) #13
colnames(output.bleed)
head(output.bleed$id) #path to image

output.modern <- read.csv("./Data/output_modern_trim.csv",
                          header = TRUE)
nrow(output.modern) #510
colnames(output.modern)
head(output.modern$id) #path to image

#keep bleed separate because magnification is different


#info about magnification for output.fossil specimen; other specimens checked in metadata_check scripts
bryo.meta <- read.csv("./Data/image_merge_txt_usingfileName_DONE_17Apr2023.csv",
                      header = TRUE,
                      sep = ";")
nrow(bryo.meta)
colnames(bryo.meta)
head(bryo.meta)

#### MANIPULATE DATA ----

## remove paths and extract just image name
output.fossil$image <- gsub("bryozoa_lab_images/",
                            "",
                            output.fossil$id)

output.fossil$image <- gsub(".jpg",
                            "",
                            output.fossil$image)
length(unique(output.fossil$image)) #1395

output.waboiv$image <- gsub("bryozoa_lab_images/WABO_IV/",
                            "",
                            output.waboiv$id)
output.waboiv$image <- gsub("/home/voje-lab/Desktop/Stegs3 WABO IV/NKLS 2-jpg/",
                            "",
                            output.waboiv$image)
output.waboiv$image <- gsub(".jpg",
                            "",
                            output.waboiv$image)
length(unique(output.waboiv$image)) #163


output.modern$image <- gsub("bryozoa_lab_images/mab-modern-jpg/",
                            "",
                            output.modern$id)
output.modern$image <- gsub(".jpg",
                            "",
                            output.modern$image)
length(unique(output.modern$image)) #76

output.bleed$image <- gsub("/home/voje-lab/Desktop/BLEED-jpg/",
                            "",
                           output.bleed$id)
output.bleed$image <- gsub(".jpg",
                            "",
                           output.bleed$image)
length(unique(output.bleed$image)) #76


meta.path.tif <- str_split(bryo.meta$path.tif, "/")
images <- c()
for(i in 1:length(meta.path.tif)){
    end <- length(meta.path.tif[[i]])
    images[i] <- meta.path.tif[[i]][end]
}
bryo.meta$old.image <- images
bryo.meta$old.image <- gsub(".tif",
                            "",
                            bryo.meta$old.image)
length(unique(bryo.meta$old.image)) #1890

##### SPECIMEN NR ----
#meta
meta.sp.nr <- bryo.meta$old.image
meta.sp.nr.sp <- str_split(meta.sp.nr, "_")
meta.sp.nr.new <- c()
for(i in 1:length(meta.sp.nr.sp)){
    meta.sp.nr.new[i] <- paste0(meta.sp.nr.sp[[i]][1], meta.sp.nr.sp[[i]][2])
}
bryo.meta$specimenNR <- meta.sp.nr.new
#check that they match
tt <- bryo.meta[,c("image", "specimenNR")]
View(tt)
length(unique(bryo.meta$specimenNR)) #905 colonies

#modern
mod.sp.nr <- output.modern$image
mod.sp.nr.sp <- str_split(mod.sp.nr, "_")
mod.sp.nr.new <- c()
for(i in 1:length(mod.sp.nr.sp)){
    mod.sp.nr.new[i] <- paste0(mod.sp.nr.sp[[i]][1], mod.sp.nr.sp[[i]][2])
}
output.modern$specimenNR <- mod.sp.nr.new
#check that they match
xx <- output.modern[,c("image", "specimenNR")]
View(xx)
length(unique(output.modern$specimenNR)) #24 colonies

#bleed
bleed.sp.nr <- output.bleed$image
bleed.sp.nr.sp <- str_split(bleed.sp.nr, "\\.")
bleed.sp.nr.new <- c()
for(i in 1:length(bleed.sp.nr.sp)){
  bleed.sp.nr.new[i] <- paste0(bleed.sp.nr.sp[[i]][3], ".", bleed.sp.nr.sp[[i]][4])
}
output.bleed$specimenNR <- bleed.sp.nr.new
#check that they match
aa <- output.bleed[,c("image", "specimenNR")]
View(aa)
length(unique(output.bleed$specimenNR)) #3 colonies

#fossil
foss.sp.nr <- output.fossil$image
foss.sp.nr.sp <- str_split(foss.sp.nr, "_")
foss.sp.nr.new <- c()
for(i in 1:length(foss.sp.nr.sp)){
    foss.sp.nr.new[i] <- paste0(foss.sp.nr.sp[[i]][1], foss.sp.nr.sp[[i]][2])
}
output.fossil$specimenNR <- foss.sp.nr.new
#check that they match
yy <- output.fossil[,c("image", "specimenNR")]
View(yy)
length(unique(output.fossil$specimenNR)) #732 colonies

#wabo iv
wabo.sp.nr <- output.waboiv$image
wabo.sp.nr.sp <- str_split(wabo.sp.nr, "_")
wabo.sp.nr.new <- c()
for(i in 1:length(wabo.sp.nr.sp)){
    wabo.sp.nr.new[i] <- paste0(wabo.sp.nr.sp[[i]][1], wabo.sp.nr.sp[[i]][2])
}
output.waboiv$specimenNR <- wabo.sp.nr.new
#check that they match
zz <- output.waboiv[,c("image", "specimenNR")]
View(zz)
length(unique(output.waboiv$specimenNR)) #50 colonies

##### FORMATION INFO -----
#modern
output.bleed$formation <- "modern"
output.bleed$form.no <- "BLEED"

output.modern$form.no <- as.numeric(str_extract(output.modern$specimenNR,
                                                "[0-9]+"))
output.modern$formation <- ""
output.modern$formation[output.modern$form.no >= 1200 &
                        output.modern$form.no <= 1299] <- "modern"

output.fossil$form.no <- as.numeric(str_extract(output.fossil$specimenNR,
                                                "[0-9]+"))
output.fossil$formation <- ""

output.fossil$formation[output.fossil$form.no <= 399] <- "NKBS"
output.fossil$formation[output.fossil$form.no >= 400 & 
                        output.fossil$form.no <= 599] <- "NKLS"
output.fossil$formation[output.fossil$form.no >= 600 & 
                        output.fossil$form.no <= 699] <- "Tewkesbury"
output.fossil$formation[output.fossil$form.no >= 700 & 
                        output.fossil$form.no <= 799] <- "SHCSBSB"
output.fossil$formation[output.fossil$form.no >= 800 & 
                        output.fossil$form.no <= 899] <- "Tainui"
output.fossil$formation[output.fossil$form.no >= 1000 & 
                        output.fossil$form.no <= 1099] <- "Upper Kai-Iwi"
output.fossil$formation[output.fossil$form.no >= 1100 & 
                        output.fossil$form.no <= 1199] <- "Tewkesbury" #formerly Waipuru
unique(output.fossil$formation)
output.fossil[output.fossil$formation == "",]
# 1200CC is Pukenni Limestone

output.waboiv$form.no <- as.numeric(str_extract(output.waboiv$specimenNR,
                                                "[0-9]+"))
output.waboiv$formation <- ""

output.waboiv$formation[output.waboiv$form.no <= 399] <- "NKBS"
output.waboiv$formation[output.waboiv$form.no >= 400 & 
                          output.waboiv$form.no <= 599] <- "NKLS"
output.waboiv$formation[output.waboiv$form.no >= 600 & 
                          output.waboiv$form.no <= 699] <- "Tewkesbury"
output.waboiv$formation[output.waboiv$form.no >= 700 & 
                          output.waboiv$form.no <= 799] <- "SHCSBSB"
output.waboiv$formation[output.waboiv$form.no >= 800 & 
                          output.waboiv$form.no <= 899] <- "Tainui"
output.waboiv$formation[output.waboiv$form.no >= 1000 & 
                          output.waboiv$form.no <= 1099] <- "Upper Kai-Iwi"
output.waboiv$formation[output.waboiv$form.no >= 1100 & 
                          output.waboiv$form.no <= 1199] <- "Tewkesbury" #formerly Waipuru
unique(output.waboiv$formation)

##### CHECK FOR DUPES -----
#fossils
output.fossil$id[1]

zz <- output.fossil[duplicated(output.fossil$box_id),]
nrow(zz) #3
zz
#882_1405_433_617
#598_1024_434_666
#880_646_517_625

output.fossil$new.id <- paste0(output.fossil$box_id, "_", output.fossil$image)
xx <- output.fossil[duplicated(output.fossil$new.id),]
nrow(xx) #0! yay!

## check for dupes
#wabo iv
output.waboiv$id[1]

zz <- output.waboiv[duplicated(output.waboiv$box_id),]
nrow(zz) #0

output.waboiv$new.id <- paste0(output.waboiv$box_id, "_", output.waboiv$image)
xx <- output.waboiv[duplicated(output.waboiv$new.id),]
nrow(xx) #0! yay!

## check for dupes
#modern
output.modern$id[1]

zz <- output.modern[duplicated(output.modern$box_id),]
nrow(zz) #0

output.modern$new.id <- paste0(output.modern$box_id, "_", output.modern$image)
xx <- output.modern[duplicated(output.modern$new.id),]
nrow(xx) #0! yay!

## check for dupes
#bleed
output.bleed$id[1]

zz <- output.bleed[duplicated(output.bleed$box_id),]
nrow(zz) #0

output.bleed$new.id <- paste0(output.bleed$box_id, "_", output.bleed$image)
xx <- output.bleed[duplicated(output.bleed$new.id),]
nrow(xx) #0! yay!


#### FILTER IMAGES ----

#modern
## remove images that are not magnifica
#images are 1200, 1203, 1205, 1216, 1220
mod.rm.img <- c("bryozoa_lab_images/mab-modern-jpg/1200_CC_1_15v_x30_BSE.jpg", #this was changed to 1222
            "bryozoa_lab_images/mab-modern-jpg/1200_CC_2_15v_x30_BSE.jpg", #this was changed to 1222
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
output.modern <- output.modern[!(output.modern$id %in% mod.rm.img),]
nrow(output.modern) #219
length(unique(output.modern$specimenNR)) #19

#remove Pukenni Limestone
output.fossil <- output.fossil[output.fossil$specimenNR != "1200CC",]
nrow(output.fossil) #6436
length(unique(output.fossil$specimenNR)) #731

#only 30 mag
unique(bryo.meta$Mag)
keep <- bryo.meta$old.image[bryo.meta$Mag == "x30"]
length(keep) #1835
length(unique(bryo.meta$specimenNR[bryo.meta$Mag == "x30"]))
setdiff(output.fossil$image, keep) #099_CV_1_15v_x35_BSE
bryo.meta$Mag[bryo.meta$image == "099_CV_1_15v_x35_BSE"] #NA, don't know why it is retained
output.fossil[output.fossil$image == "099_CV_1_15v_x35_BSE",] #5 
output.fossil <- output.fossil[output.fossil$image %in% keep,]
nrow(output.fossil) #6431; lost 5 images, all 099_CV_1_15v_x35_BSE

unique(df.filter$Mag) #only x30
filter.path.tif <- str_split(fossil.filter$path.tif, "/")
images.filter <- c()
for(i in 1:length(filter.path.tif)){
  end <- length(filter.path.tif[[i]])
  images.filter[i] <- filter.path.tif[[i]][end]
}
fossil.filter$old.image <- images.filter
fossil.filter$old.image <- gsub(".tif",
                            "",
                            fossil.filter$old.image)
setdiff(fossil.filter$old.image, output.fossil$image)
setdiff(output.fossil$image, fossil.filter$old.image) #no diff
nrow(output.fossil)
output.fossil <- output.fossil[output.fossil$image %in% fossil.filter$old.image,]
nrow(output.fossil) #no diff, as expected

##remove other:
# 005_CV_1: have two of them (NKBS)
# 099_CV_1: mag is 35 (NKBS)
foss.rm.img <- c("005_CV_1_15v_x30", "099_CV_1_15v_x35_BSE")
output.fossil <- output.fossil[!(output.fossil$image %in% foss.rm.img),]

## filter bleed
output.bleed <- output.bleed[output.bleed$image %in% bleed.filter$image,]
nrow(output.bleed) #11, lost 2 colonies

#### COMBINE DATASETS AND WRITE OUT DATASET ----
#make columns match
setdiff(colnames(output.fossil), colnames(output.modern))
setdiff(colnames(output.modern), colnames(output.fossil))

setdiff(colnames(output.fossil), colnames(output.waboiv))
setdiff(colnames(output.waboiv), colnames(output.fossil))

setdiff(colnames(output.fossil), colnames(output.bleed))
setdiff(colnames(output.bleed), colnames(output.fossil))

images.df <- as.data.frame(rbind(output.fossil,
                                 output.waboiv,
                                 output.modern,
                                 output.bleed))
nrow(images.df) #7254

write.csv(images.df,
          "./Data/images.filtered_27May2024.csv",
          row.names = FALSE)

# #filter images
# colnames(df.filter)
# head(df.filter$fileName) #new file name
# head(df.filter$path.tif)
# #extract only the image name, not the entire path; helps match df.filter and images.meta
# #merge by images.meta$image
# #path.tif contains old file name
# df.filter$fileName.old <- c()
# for(i in 1:nrow(df.filter)){
#     df.filter$fileName.old[i] <- str_split(df.filter$path.tif[i], "/")[[1]][length(str_split(df.filter$path.tif[i], "/")[[1]])]
# }
# df.filter$fileName.old <- gsub(".tif",
#                                "",
#                                df.filter$fileName.old)
# 
# setdiff(meta.images.fossil$image, df.filter$fileName.old)
# #only one different is: 099_CV_1_15v_x35_BSE
# nrow(meta.images.fossil) #6443
# images.filter.fossil <- meta.images.fossil[meta.images.fossil$image %in% df.filter$fileName.old,]
# nrow(images.filter.fossil) #6438 #lost 5, likely all 099_CV_1_15v_x35_BSE
# length(unique(images.filter.fossil$fileName.tif)) #1394
# 



#check if 5 of each colony
path.trim <- gsub(wabo.iv$id,
                  pattern = "bryozoa_lab_images/WABO_IV/",
                  replacement = "")
path.split <- str_split(path.trim,
                        pattern = "_")

mag <- c()
col.id <- c()
col <- c()
side <- c()
for(i in 1:length(path.split)){
  mag[i] <- path.split[[i]][5]
  col[i] <- path.split[[i]][1]
  side[i] <- path.split[[i]][2]
  col.id[i] <- paste0(path.split[[i]][1], path.split[[i]][2])
}
mag == "x30" #all true

colonies <- as.data.frame(cbind(col, side, col.id))
sort(table(colonies$col.id))
few <- c("076iCV", "074iCV", "1028CV", "1029CV", "843F",
         "065iCC", "073iaCV", "073ibCV")

colonies$form.no <- as.numeric(str_extract(colonies$col,
                                           "[0-9]+"))
colonies$formation <- ""

colonies$formation[colonies$form.no <= 399] <- "NKBS"
colonies$formation[colonies$form.no >= 400 & 
                     colonies$form.no <= 599] <- "NKLS"
colonies$formation[colonies$form.no >= 600 & 
                     colonies$form.no <= 699] <- "Tewkesbury"
colonies$formation[colonies$form.no >= 700 & 
                     colonies$form.no <= 799] <- "SHCSBSB"
colonies$formation[colonies$form.no >= 800 & 
                     colonies$form.no <= 899] <- "Tainui"
colonies$formation[colonies$form.no >= 1000 & 
                     colonies$form.no <= 1099] <- "Upper Kai-Iwi"
colonies$formation[colonies$form.no >= 1100 & 
                     colonies$form.no <= 1199] <- "Tewkesbury" #formerly Waipuru
unique(colonies$formation)
colonies[colonies$formation == "",]

colonies[!(colonies$col.id %in% few),] %>%
  group_by(formation) %>%
  summarise(n = length(unique(col.id)))

#modern
output.modern$id[1]

aa <- output.modern[duplicated(output.modern$box_id),]
nrow(aa) #0

output.modern$new.id <- paste0(output.modern$box_id, "_", output.modern$image)
vv <- output.modern[duplicated(output.modern$new.id),]
nrow(vv) #0! yay!
