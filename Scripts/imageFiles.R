## Meghan A. Balk
## meghan.balk@nhm.uio.no

## This code:
## 1) extracts file names from the zip file of Steginoporella images
## 2) creates a csv file with the information parsed

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
list = list.files(path = "/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images",
                  full.names = TRUE,
                  recursive = TRUE)

#get rid of abstract book, I think
list.rm <- list[!grepl("/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images/Stegs/ISME 14 ABSTRACT BOOK", list)]

#### CREATE CSV ----

listPath <- unlist(list.rm)
length(listPath) #3779

##### PARSE FILE NAMES -----

list.trim <- gsub(list.rm,
                  pattern = "/Users/mab/Library/CloudStorage/Dropbox/Rocks-Paradox/Bryozoans/Stegino images/",
                  replacement = "")

list.parse <- str_split(list.trim,
                        pattern = "/")

#first folder is either Sara (folder name Sara) or Mali (folder names Stegs and Stegs2)
#subfolder folder, when given, is the formation or grouping

folder <- c()
formation <- c()
fileName <- c()
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
    fileName[i] <- list.parse[[i]][2]
  }
  else if(isTRUE(endsWith(list.parse[[i]][2], ".tif"))){
    fileName[i] <- list.parse[[i]][2]
  }
  else{
    fileName[i] <- list.parse[[i]][3]
  }
  if(isTRUE(endsWith(fileName[i], ".txt"))){
    ext[i] <- "txt"
  }
  else{
    ext[i] <- "tif"
  }
}

##### PARSE IMAGE NAME -----

image <- str_extract(fileName, pattern = "[^.]+")

ImageName.list <- str_split(image,
                            pattern = "_")

specimenNR <- c()
number <- c()
colonyCurve <- c()
pictureNumber <- c()
AV <- c()
magnification <- c()
backscatter <- c()

for(i in 1:length(ImageName.list)){
  number[i] <- ImageName.list[[i]][1]
  colonyCurve[i] <- ImageName.list[[i]][2]
  pictureNumber[i] <- ImageName.list[[i]][3]
  AV[i] <- ImageName.list[[i]][4]
  magnification[i] <- str_extract(ImageName.list[[i]][5],
                                  pattern = "[^.]+")
  backscatter[i] <- str_extract(ImageName.list[[i]][6],
                                pattern = "[^.]+")
  specimenNR[i] <- paste0(ImageName.list[[i]][1], ImageName.list[[i]][2])
}

##### COMBINE & WRITE CSV ----

df.list <- data.frame(path = listPath,
                      folder = folder,
                      formation = formation,
                      image = image,
                      ext = ext,
                      fileName = fileName,
                      specimenNR = specimenNR,
                      number = number,
                      colonyCurve = colonyCurve,
                      pictureNumber = pictureNumber,
                      AV = AV,
                      magnification = magnification,
                      backscatter = backscatter,
                      stringsAsFactors = FALSE)

nrow(df.list) #3779
nrow(df.list[df.list$ext == "tif",]) #1890, an extra image

df.list$newFormation <- ""
df.list$newFormation[grepl("^0", df.list$specimenNR)] <- "NKBS"
df.list$newFormation[grepl("^1", df.list$specimenNR)] <- "NKBS"
df.list$newFormation[grepl("^2", df.list$specimenNR)] <- "NKBS"
df.list$newFormation[grepl("^3", df.list$specimenNR)] <- "NKBS"
df.list$newFormation[grepl("^4", df.list$specimenNR)] <- "NKLS"
df.list$newFormation[grepl("^5", df.list$specimenNR)] <- "NKLS"
df.list$newFormation[grepl("^6", df.list$specimenNR)] <- "Tewkesbury"
df.list$newFormation[grepl("^7", df.list$specimenNR)] <- "SHCSBSB"
df.list$newFormation[grepl("^8", df.list$specimenNR)] <- "Tainui"
df.list$newFormation[grepl("^10", df.list$specimenNR)] <- "Upper Kai-Iwi"
df.list$newFormation[grepl("^11", df.list$specimenNR)] <- "Waipuru"

#write.csv(df.list,
#          "./Data/imageList.csv",
#          row.names = FALSE)

