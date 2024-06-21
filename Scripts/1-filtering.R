
# Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
# 1) adds in formation info
# 2) filters images
# 3) makes sure that all datasets have:
# formation, specimenNR, sources (i.e., BLEED/EPA, NIWA, pdt), and magnification

##Notes
# use all magnification
# for BLEED images, use lower magnification as higher magnification are zoomed-
# in for the same part of the colony (see column "use" and select "Y")
#

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

##### FILTERING TABLES -----
fossil.filter <- read.csv("Data/filteredImages.csv", sep = ";", header = TRUE)
#bleed.filter <- read.csv("Data/BLEED.image.filter.csv", header = TRUE)
#keeping all wabo iv, so no filter file needed
#manually filtering modern, so no filter file needed

##### IMAGE METADATA -----
#info about magnification for output.fossil specimen; other specimens checked in metadata_check scripts
bryo.meta <- read.csv("./Data/metadata/image_merge_txt_usingfileName_DONE_17Apr2023.csv",
                      header = TRUE,
                      sep = ";")
nrow(bryo.meta)
colnames(bryo.meta)
head(bryo.meta)

klv.meta <- read.csv("Data/metadata/KLV_metadata.csv",
                     header = TRUE)

bleed.meta <- read.csv("Data/metadata/BLEED_metadata.csv",
                       header = TRUE)

niwa.meta <- read.csv("Data/metadata/NIWA_metadata.csv",
                      header = TRUE)

##### OUTPUTS -----
## FOSSILS
output.fossil <- read.csv("./Data/steginator_magnifica_outputs/output_4Aug2023_done.csv", header = TRUE)
nrow(output.fossil) #6443
colnames(output.fossil)
head(output.fossil$id) #path to image

output.waboiv_1 <- read.csv("Data/steginator_magnifica_outputs/output_waboiv.csv", header = TRUE)
output.waboiv_2 <- read.csv("Data/steginator_magnifica_outputs/output.nkls.2.csv", header = TRUE)

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

## MODERN

output.bleed1 <- read.csv("Data/steginator_magnifica_outputs/output.bleed.csv", header = TRUE)

nrow(output.bleed1) #11
colnames(output.bleed1)
head(output.bleed1$id) #path to image

output.bleed2 <- read.csv("Data/steginator_magnifica_outputs/output_bleed2.csv", header = TRUE)

nrow(output.bleed2) #11
colnames(output.bleed2)
head(output.bleed2$id) #path to image

#### COMBINE BLEED
output.bleed <- as.data.frame(rbind(output.bleed1, output.bleed2))

##remove duplicate IDs
duplicated(output.bleed$box_id)
output.bleed <- output.bleed[!duplicated(output.bleed$box_id),]

nrow(output.bleed) #17
colnames(output.bleed)
head(output.bleed$id) #path to image

output.modern <- read.csv("./Data/steginator_magnifica_outputs/output_modern_trim.csv",
                          header = TRUE)
nrow(output.modern) #510
colnames(output.modern)
head(output.modern$id) #path to image

output.niwa <- read.csv("Data/steginator_magnifica_outputs/output_niwa2.csv",
                        header = TRUE)
nrow(output.niwa) #197
colnames(output.niwa)
head(output.niwa$id) #path to image

output.klv <- read.csv("Data/steginator_magnifica_outputs/output_klv.csv",
                       header = TRUE)
nrow(output.klv) #7
colnames(output.klv)
head(output.klv$id) #path to image

#keep bleed separate because magnification is different

#### MANIPULATE DATA ----

##### IMAGE/FILE NAME ----
## remove paths and extract just image name
###### FOSSIL ------
output.fossil$image <- gsub("bryozoa_lab_images/",
                            "",
                            output.fossil$id)

output.fossil$image <- gsub(".jpg",
                            "",
                            output.fossil$image)
head(output.fossil$image)
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
head(output.waboiv$image)
length(unique(output.waboiv$image)) #163

###### MODERN ------
output.modern$image <- gsub("bryozoa_lab_images/mab-modern-jpg/",
                            "",
                            output.modern$id)
output.modern$image <- gsub(".jpg",
                            "",
                            output.modern$image)
head(output.modern$image)
length(unique(output.modern$image)) #76

output.bleed$image <- gsub("/home/voje-lab/Desktop/BLEED-jpg/",
                            "",
                           output.bleed$id)
output.bleed$image <- gsub(".jpg",
                            "",
                           output.bleed$image)
head(output.bleed$image)
length(unique(output.bleed$image)) #10

output.klv$image <- gsub("/home/voje-lab/Desktop/KLV_stegs-jpg/",
                         "",
                         output.klv$id)
output.klv$image <- gsub(".jpg",
                         "",
                         output.klv$image)
head(output.klv$image)
length(unique(output.klv$image)) #3

output.niwa$image <- gsub("/home/voje-lab/Desktop/NIWA_magnifica_SEM-jpg/Steginoporella magnifica/",
                         "",
                         output.niwa$id)
output.niwa$image <- gsub(".jpg",
                         "",
                         output.niwa$image)
head(output.niwa$image)
length(unique(output.niwa$image)) #3

###### METADATA ------
#only an issue for bryo.meta
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

##### ADD SOURCE -----
output.waboiv$source <- "BLEED/EPA"
output.fossil$source <- "BLEED/EPA"
output.modern$source <- "BLEED/EPA"
output.bleed$source <- "BLEED/EPA"
output.niwa$source <- "NIWA"
output.klv$source <- c("pdt", "pdt", "pdt", "pdt", "pdt", "NIWA", "NIWA")

##### SPECIMEN NR ----
###### METADATA -----
#only for bryo.meta
meta.sp.nr <- bryo.meta$old.image
meta.sp.nr.sp <- str_split(meta.sp.nr, "_")
meta.sp.nr.new <- c()
meta.image.no <- c()
for(i in 1:length(meta.sp.nr.sp)){
    meta.sp.nr.new[i] <- paste0(meta.sp.nr.sp[[i]][1], meta.sp.nr.sp[[i]][2])
    meta.image.no[i] <- meta.sp.nr.sp[[i]][3]
}
bryo.meta$specimenNR <- meta.sp.nr.new
bryo.meta$image.no <- meta.image.no
#check that they match
View(bryo.meta[,c("image", "specimenNR", "image.no")])
length(unique(bryo.meta$specimenNR)) #905 colonies

###### FOSSIL ------
foss.sp.nr <- output.fossil$image
foss.sp.nr.sp <- str_split(foss.sp.nr, "_")
foss.sp.nr.new <- c()
foss.image.no <- c()
for(i in 1:length(foss.sp.nr.sp)){
    foss.sp.nr.new[i] <- paste0(foss.sp.nr.sp[[i]][1], foss.sp.nr.sp[[i]][2])
    foss.image.no[i] <- foss.sp.nr.sp[[i]][3]
}
output.fossil$specimenNR <- foss.sp.nr.new
output.fossil$image.no <- foss.image.no
#check that they match
View(output.fossil[,c("image", "specimenNR", "image.no")])
length(unique(output.fossil$specimenNR)) #732 colonies

#wabo iv
wabo.sp.nr <- output.waboiv$image
wabo.sp.nr.sp <- str_split(wabo.sp.nr, "_")
wabo.sp.nr.new <- c()
wabo.image.no <- c()
for(i in 1:length(wabo.sp.nr.sp)){
    wabo.sp.nr.new[i] <- paste0(wabo.sp.nr.sp[[i]][1], wabo.sp.nr.sp[[i]][2])
    wabo.image.no[i] <- wabo.sp.nr.sp[[i]][3]
}
output.waboiv$specimenNR <- wabo.sp.nr.new
output.waboiv$image.no <- wabo.image.no
#check that they match
View(output.waboiv[,c("image", "specimenNR", "image.no")])
length(unique(output.waboiv$specimenNR)) #50 colonies

###### MODERN ------
mod.sp.nr <- output.modern$image
mod.sp.nr.sp <- str_split(mod.sp.nr, "_")
mod.sp.nr.new <- c()
mod.image.no <- c()
for(i in 1:length(mod.sp.nr.sp)){
    mod.sp.nr.new[i] <- paste0(mod.sp.nr.sp[[i]][1], mod.sp.nr.sp[[i]][2])
    mod.image.no[i] <- mod.sp.nr.sp[[i]][3]
}
output.modern$specimenNR <- mod.sp.nr.new
output.modern$image.no <- mod.image.no
#check that they match
View(output.modern[,c("image", "specimenNR", "image.no")])
length(unique(output.modern$specimenNR)) #24 colonies

#bleed
bleed.sp.nr <- output.bleed$image
bleed.sp.nr.sp <- str_split(bleed.sp.nr, "\\.")
bleed.sp.nr.new <- c()
bleed.image.no <- c()
for(i in 1:length(bleed.sp.nr.sp)){
  bleed.sp.nr.new[i] <- bleed.sp.nr.sp[[i]][4]
  bleed.image.no[i] <- paste0(bleed.sp.nr.sp[[i]][1], ".", bleed.sp.nr.sp[[i]][2])
}
output.bleed$specimenNR <- bleed.sp.nr.new
output.bleed$image.no <- bleed.image.no
#check that they match
View(output.bleed[,c("image", "specimenNR", "image.no")])
length(unique(output.bleed$specimenNR)) #7 colonies

#niwa
##there are multiple colonies on a specimen, so here, specimen nr is really colony number
niwa.sp.nr <- output.niwa$image
niwa.sp.nr.sp <- str_split(niwa.sp.nr, " ")
niwa.sp.nr.new <- c()
niwa.image.no <- c()
for(i in 1:length(niwa.sp.nr.sp)){
    niwa.sp.nr.new[i] <- paste0(niwa.sp.nr.sp[[i]][3], " ", niwa.sp.nr.sp[[i]][4], 
                                " ", str_extract(niwa.sp.nr.sp[[i]][6], "[a-zA-Z]"))
    niwa.image.no[i] <- str_extract(niwa.sp.nr.sp[[i]][6], "[0-9]+")
} #has NA but too lazy to figure it out
niwa.sp.nr.new <- gsub(" NA", replacement = "", niwa.sp.nr.new)
output.niwa$specimenNR <- niwa.sp.nr.new
output.niwa$image.no <- niwa.image.no
#check that they match
View(output.niwa[,c("image", "specimenNR", "image.no")])
length(unique(output.niwa$specimenNR)) #11 colonies
##need to fix one image
output.niwa$specimenNR[47:49] <- "NIWA 122574"
output.niwa$image.no[47:49] <- 2

#klv
output.klv$image
## do by hand, because not worth the effort
output.klv$specimenNR <- c("pdt20893", "pdt20893",
                            "pdt21191", "pdt21191", "pdt21191",
                            "NIWA98542", "NIWA98542")
output.klv$image.no <- rep("1", nrow(output.klv))
View(output.klv[,c("image", "specimenNR", "image.no")])
length(unique(output.klv$specimenNR)) #3 colonies

##### ADD IN METADATA -----
###### FOSSIL ------
#bryo.meta does not include wabo iv
#i have manually checked the metadata files for wabo iv and all mag is 30 and correct
bryo.meta.trim <- bryo.meta %>%
    select("image", "Mag", "specimenNR", "image.no")
nrow(bryo.meta)
nrow(output.fossil)
fossil <- merge(output.fossil,
                bryo.meta.trim,
                by = c("image", "specimenNR", "image.no"),
                all.x = TRUE, all.y = FALSE)
nrow(fossil)
#colnames(fossil)
#View(fossil[, c(1, 2, 3, 57, 58)])

##waboiv
wabo <- output.waboiv
wabo$Mag <- rep("x30", nrow(wabo)) #already checked

###### MODERN ------
##modern
mod <- output.modern 
mod$Mag <- rep("x30", nrow(mod)) #should all have mag 30

##bleed
bleed.meta.trim <- bleed.meta %>%
    select("file_name", "SEM_no", "BLEED_code", "magnification", "use")
colnames(bleed.meta.trim) <- c("image", "image.no", "specimenNR", "Mag", "use")
nrow(bleed.meta.trim)
nrow(output.bleed)
bleed <- merge(output.bleed,
               bleed.meta.trim,
               by = c("image", "specimenNR", "image.no"),
               all.x = TRUE, all.y = FALSE)
nrow(bleed)
##need to check and fix; should be 294 not 293
bleed$Mag[bleed$image.no == "MHR.293"] <- "x80"
bleed$use[bleed$image.no == "MHR.293"] <- "N"

bleed <- bleed[bleed$use == "Y",]
bleed <- bleed[, -59]

##niwa
niwa.meta.trim <- niwa.meta %>%
    select("File_name", "specimen_no", "image_no", "mag", "use")
colnames(niwa.meta.trim) <- c("image", "specimenNR", "image.no", "Mag", "use")
nrow(niwa.meta.trim)
nrow(output.niwa)
niwa <- merge(output.niwa,
              niwa.meta.trim,
              by = c("image", "specimenNR", "image.no"),
              all.x = TRUE, all.y = FALSE)
nrow(niwa)
View(niwa[, c(1, 2, 3, 57, 58, 59)])

niwa <- niwa[niwa$use == "Y",]
niwa <- niwa[, -59]

##klv
klv.meta.trim <- klv.meta %>%
    select("File_name", "specimen_no", "image_no", "mag")
#actually going to do by hand because it's a mixed dataset
klv <- output.klv
klv$Mag <- c("x40", "x40", "x50", "x50", "x50", "x50", "x50")

##### SCALE -----
# all lab photos have the following scales:
# x30 = 0.606
# x35 = 0.705
# x40 = 0.825
# x50 = 0.908
# x60 = 1.089
# x80 = 1.618
# all niwa photos have the following scales:
# x30 = 0.242
# x50 = 0.404
# x60 = 0.484
# all pdt photos have the following scales:
# x40 = 0.302
# x50 = 0.376

###### FOSSIL ------
unique(fossil$Mag)
unique(fossil$image[is.na(fossil$Mag)])
#need to fix a couple; checked and all are x30 except 005iF image 1
fossil$Mag[fossil$image %in% unique(fossil$image[is.na(fossil$Mag)])] <- "x30"
fossil$scale <- 0.606

unique(wabo$Mag)
wabo$scale <- 0.606

###### MODERN ------
unique(mod$Mag)
mod$scale <- 0.606

unique(bleed$Mag)
bleed$scale <- ""
bleed$scale[bleed$Mag == "x50"] <- 0.908
bleed$scale[bleed$Mag == "x40"] <- 0.825

unique(niwa$Mag)
niwa$scale <- ""
niwa$scale[niwa$Mag == "x30"] <- 0.242
niwa$scale[niwa$Mag == "x50"] <- 0.404
niwa$scale[niwa$Mag == "x60"] <- 0.484

unique(klv$Mag)
#again by hand
klv$scale <- c(0.302, 0.302, 0.376, 0.376, 0.376,  0.404,  0.404)

##### LOCALITY / TIME INFO -----
###### FOSSIL ------
fossil$form.no <- as.numeric(str_extract(fossil$specimenNR,
                                         "[0-9]+"))
fossil$formation <- ""

fossil$formation[fossil$form.no <= 399] <- "NKBS"
fossil$formation[fossil$form.no >= 400 & 
                     fossil$form.no <= 599] <- "NKLS"
fossil$formation[fossil$form.no >= 600 & 
                     fossil$form.no <= 699] <- "Tewkesbury"
fossil$formation[fossil$form.no >= 700 & 
                     fossil$form.no <= 799] <- "SHCSBSB"
fossil$formation[fossil$form.no >= 800 & 
                     fossil$form.no <= 899] <- "Tainui"
fossil$formation[fossil$form.no >= 1000 & 
                     fossil$form.no <= 1099] <- "Upper Kai-Iwi"
fossil$formation[fossil$form.no >= 1100 & 
                     fossil$form.no <= 1199] <- "Tewkesbury" #formerly Waipuru
unique(fossil$formation)
fossil[fossil$formation == "",]
# 1200CC is Pukenni Limestone
fossil <- fossil[fossil$specimenNR != "1200CC",] #removing it

wabo$form.no <- as.numeric(str_extract(wabo$specimenNR,
                                       "[0-9]+"))
wabo$formation <- ""

wabo$formation[wabo$form.no <= 399] <- "NKBS"
wabo$formation[wabo$form.no >= 400 & 
                   wabo$form.no <= 599] <- "NKLS"
wabo$formation[wabo$form.no >= 600 & 
                   wabo$form.no <= 699] <- "Tewkesbury"
wabo$formation[wabo$form.no >= 700 & 
                   wabo$form.no <= 799] <- "SHCSBSB"
wabo$formation[wabo$form.no >= 800 & 
                   wabo$form.no <= 899] <- "Tainui"
wabo$formation[wabo$form.no >= 1000 & 
                   wabo$form.no <= 1099] <- "Upper Kai-Iwi"
wabo$formation[wabo$form.no >= 1100 & 
                   wabo$form.no <= 1199] <- "Tewkesbury" #formerly Waipuru
unique(wabo$formation)

###### MODERN ------
mod$formation <- "modern"
mod$form.no <- as.numeric(str_extract(mod$specimenNR,
                                      "[0-9]+"))

bleed$formation <- "modern"
bleed$form.no <- "BLEED"

niwa$formation <- "modern"
niwa$form.no <- "NIWA"

#again, by hand
klv$formation <- c("Tainui", "Tainui", "modern", "modern", "modern", "modern", "modern")
klv$form.no <- "klv"

##### CHECK FOR DUPES -----
###### FOSSILS ------
fossil$id[1]

zz <- fossil[duplicated(fossil$box_id),]
nrow(zz) #3
fossil[fossil$box_id %in% zz$box_id,]
#882_1405_433_617
#598_1024_434_666
#880_646_517_625
#has same box_id but different images

fossil$new.id <- paste0(fossil$box_id, "_", fossil$image)
xx <- fossil[duplicated(fossil$new.id),]
nrow(xx) #0! yay!

## check for dupes
#wabo iv
wabo$id[1]
nrow(wabo[duplicated(wabo$box_id),]) #0

wabo$new.id <- paste0(wabo$box_id, "_", wabo$image)
nrow(wabo[duplicated(wabo$new.id),])

###### MODERN ------
nrow(mod[duplicated(mod$box_id),]) #0

mod$new.id <- paste0(mod$box_id, "_", mod$image)
nrow(mod[duplicated(mod$new.id),])

#bleed
nrow(bleed[duplicated(bleed$box_id),])

bleed$new.id <- paste0(bleed$box_id, "_", bleed$image)
nrow(bleed[duplicated(bleed$new.id),])

#niwa
nrow(niwa[duplicated(niwa$box_id),])

niwa$new.id <- paste0(niwa$box_id, "_", niwa$image)
nrow(niwa[duplicated(niwa$new.id),])

#klv
nrow(klv[duplicated(klv$box_id),])

klv$new.id <- paste0(klv$box_id, "_", klv$image)
nrow(klv[duplicated(klv$new.id),])

#### FILTER IMAGES ----

##### FOSSIL -----
# remove 005_CV_1_15v_x30
fossil <- fossil[fossil$image != "005_CV_1_15v_x30",]


unique(fossil.filter$Mag) #only x30
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
setdiff(fossil.filter$old.image, fossil$image) #expected as have waboiv in here now
setdiff(fossil$image, fossil.filter$old.image) #099_CV_1_15v_x35_BSE, that's fine

##### MODERN -----
#already removed some based on the column "use" above

## remove images that are not magnifica
#images are 1200, 1203, 1205, 1216, 1220
mod.rm.img <- c("1200CC", #this was changed to 1222
                "1203CC",
                "1205CC",
                "1216CC",
                "1220CC")
mod <- mod[!(mod$specimenNR %in% mod.rm.img),]
nrow(mod) #219
length(unique(mod$specimenNR)) #19

#### COMBINE DATASETS AND WRITE OUT DATASET ----
#make columns match
setdiff(colnames(fossil), colnames(mod))
setdiff(colnames(mod), colnames(fossil))

setdiff(colnames(fossil), colnames(wabo))
setdiff(colnames(wabo), colnames(fossil))

setdiff(colnames(fossil), colnames(bleed))
setdiff(colnames(bleed), colnames(fossil))

setdiff(colnames(fossil), colnames(niwa))
setdiff(colnames(niwa), colnames(fossil))

setdiff(colnames(fossil), colnames(klv))
setdiff(colnames(klv), colnames(fossil))

images.df <- as.data.frame(rbind(fossil,
                                 wabo,
                                 mod,
                                 bleed,
                                 niwa,
                                 klv))
nrow(images.df) #7391

write.csv(images.df,
          "./Data/images.filtered_21Jun2024.csv",
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

#### MAKE A SEPARATE FILE FOR MODERN AND FOSSIL OF SAMPLING INFO ----
# do by hand once have the final list of images