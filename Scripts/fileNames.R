## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
## 1) extracts file names from the zip file of Steginoporella images
## 2) creates a csv file with the information parsed
## 3) compares the files to the bryozoan metadata file

#### LOAD PACKAGES ----
require(stringr)
require(dplyr)

#### EXTRACT FILE NAMES ----

##https://stackoverflow.com/questions/54510134/getting-list-of-file-names-in-a-directory
# To list all files of a folder in a list variable including files 
# from sub-folders. The code below gets the full path of files not just names.
#list = list.files(path = full_path_to_directory ,full.names=TRUE,recursive=TRUE)
# To get names of all files from their corresponding paths in all_names variable.
#all_names = basename(list)
# To write all_names variable to a CSV file.
#write.csv(all_names, "test.csv")

## get folder names
list = list.files(path = "/Users/meghanabalk/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images",
                  full.names = TRUE,
                  recursive = TRUE)
list.trim <- gsub(list,
                  pattern = "/Users/meghanabalk/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images/",
                  replacement = "")

#get rid of abstract book, I think
list.rm <- list.trim[!grepl("Stegs/ISME 14 ABSTRACT BOOK", list.trim)]

#### CREATE CSV ----

##### PARSE FILE NAMES -----
list.parse <- str_split(list.rm,
                        pattern = "/")

#first folder is either Sara (folder name Sara) or Mali (folder names Stegs and Stegs2)
#subfolder folder, when given, is the formation or grouping

folder <- c()
formation <- c()
image <- c()
ext <- c()

for(i in 1:length(list.parse)){
  folder[i] <- list.parse[[i]][1]
  if(isTRUE(endsWith(list.parse[[i]][2], ".txt"))){
    formation[i] <- "NONE"
  }
  else if(isTRUE(endsWith(list.parse[[i]][2], ".tif"))){
    formation[i] <- "NONE"
  }
  else{
    formation[i] <- list.parse[[i]][2]
  }
  if(isTRUE(endsWith(list.parse[[i]][2], ".txt"))){
    image[i] <- list.parse[[i]][2]
  }
  else if(isTRUE(endsWith(list.parse[[i]][2], ".tif"))){
    image[i] <- list.parse[[i]][2]
  }
  else{
    image[i] <- list.parse[[i]][3]
  }
  if(isTRUE(endsWith(image[i], ".txt"))){
    ext[i] <- "txt"
  }
  else{
    ext[i] <- "tif"
  }
}

##### PARSE IMAGE NAME -----

imageName <- paste(str_extract(image, pattern = "[^.]+"), formation, sep = "_")

imageName.list <- str_split(image,
                            pattern = "_")

specimenNR <- c()
number <- c()
colonyCurve <- c()
pictureNumber <- c()
AV <- c()
magnification <- c()
backscatter <- c()

for(i in 1:length(imageName.list)){
  number[i] <- imageName.list[[i]][1]
  colonyCurve[i] <- imageName.list[[i]][2]
  pictureNumber[i] <- imageName.list[[i]][3]
  AV[i] <- imageName.list[[i]][4]
  magnification[i] <- str_extract(imageName.list[[i]][5],
                                  pattern = "[^.]+")
  backscatter[i] <- str_extract(imageName.list[[i]][6],
                                pattern = "[^.]+")
  specimenNR[i] <- paste0(imageName.list[[i]][1], imageName.list[[i]][2])
}

##### COMBINE & WRITE CSV ----

df.list <- data.frame(folder = folder,
                      formation = formation,
                      image = image,
                      ext = ext,
                      imageName = imageName,
                      specimenNR = specimenNR,
                      number = number,
                      colonyCurve = colonyCurve,
                      pictureNumber = pictureNumber,
                      AV = AV,
                      magnification = magnification,
                      backscatter = backscatter,
                      stringsAsFactors = FALSE)

nrow(df.list) #3779
nrow(df.list[df.list$ext == "tif",])

write.csv(df.list,
          "imageList.csv",
          row.names = FALSE)

#### COMPARE TO METADATA FILE ----

##### LOAD DATA -----
bryo.meta <- read.csv("Imaged Steginoporella magnifica specimens.csv", header = TRUE)
bryo.meta$SPECIMEN.NR <- gsub(bryo.meta$SPECIMEN.NR,
                              pattern = " ", 
                              replacement = "")
nrow(bryo.meta) #880
nrow(bryo.meta[!duplicated(bryo.meta$SPECIMEN.NR),]) #777

##### COMPARE TOTALS -----

tots <- df.list %>%
  group_by(specimenNR, ext) %>%
  summarise(N = n()) %>%
  as.data.frame()

nrow(tots) #1810
#1033 more images

length(setdiff(tots$specimenNR, bryo.meta$SPECIMEN.NR)) #136 in df.list that are not in bryo.meta
length(setdiff(bryo.meta$SPECIMEN.NR, tots$specimenNR)) #8 in bryo.meta that are not in df.list

## look at 8 in bryo.meta that are not in df.list
setdiff(bryo.meta$SPECIMEN.NR, tots$specimenNR) 
#e.g., bryo.meta has 204CV, df.list has 204CC
#232CV in bryo.meta is 232CC in df.list
#005CC in bryo.meta is 005CV in df.list
#0367iCV does not exist in df.list
#435CC does not exist in df.list
#466CV in bryo.meta is 466CC in df.list
#568CC in bryo.meta is 568CV in df.list
#717iCV does not exist in df.list, but 717CV does

#NR OF PICS in bryo.meta don't match the total number of images taken
max(bryo.meta$NR.OF.PICS) #4
max(tots$N) #5
tots$specimenNR[which.max(tots$N)] #237CC
df.list[df.list$specimenNR == "237CC",]
#should be only 3 images
#has 237_CC_2_15v_x40_BSE replicated
#has 237_CC_1_15v_x30_BSE_BTK2

## look at 136 in df.list that are not in bryo.meta
setdiff(tots$specimenNR, bryo.meta$SPECIMEN.NR)


extra.df <- setdiff(tots$specimenNR, bryo.meta$SPECIMEN.NR)
extra.bryo <- setdiff(bryo.meta$SPECIMEN.NR, tots$specimenNR)

df.trim <- df.list[df.list$specimenNR %in% extra.df &
                   df.list$ext == "tif",]
df.trimmed <- df.trim[!duplicated(df.trim[c('specimenNR')]), ]

dates <- c(rep("NA", length(extra.df)), bryo.meta$DATE[bryo.meta$SPECIMEN.NR %in% extra.bryo])
folder <- c(df.trimmed$folder, rep("NA", length(extra.bryo)))
formation <- c(df.trimmed$formation, rep("NA", length(extra.bryo)))
fromDataset <- c(rep("imageFiles", length(extra.df)), rep("metadata", length(extra.bryo)))
imageDiff <- c(extra.df, extra.bryo)
errors.df <- cbind(fromDataset, imageDiff, folder, formation, dates)

write.csv(errors.df,
          "errors.csv",
          row.names = FALSE)

##### COMPARE FORMATION NAMES -----
#don't have the same formations
bryo.forms <- unique(bryo.meta$FORMATION)
image.forms <- unique(df.list$formation)
#non-standardized names
inDataset <- c(rep("imageFiles", length(image.forms)), rep("metadata", length(bryo.forms)))
formations <- c(image.forms, bryo.forms)

formations.df <- cbind(inDataset, formations)

write.csv(formations.df,
          "formations.csv",
          row.names = FALSE)

## NEED META DATA FOR FORMATION
#include locality, age range
